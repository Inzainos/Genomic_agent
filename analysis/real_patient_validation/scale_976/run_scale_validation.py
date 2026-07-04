"""Scale validation runner (976 TCGA-BRCA patients) — as run 2026-07-03.
Record of the batched GDC POST /data pipeline. Paths (WORK, DB) point at the
run scratchpad; adjust them to re-run. Writes scale_results.jsonl (per patient)
and scale_summary.json (aggregate). See SCALE_VALIDATION.md.
"""
import json, os, sys, io, tarfile, time, sqlite3, urllib.request, statistics
os.environ["SNT_LOG_LEVEL"]="ERROR"
DB="/tmp/snt_scale.db"; os.environ["SNT_DB_PATH"]=DB
GA="/home/user/The-shadow-Node-Theory/genomic_agent"
sys.path.insert(0,GA+"/genomic_database"); sys.path.insert(0,GA+"/agent_core")

WORK="/tmp/claude-0/-home-user-The-shadow-Node-Theory/f6139116-5ed3-5a9c-8cc9-f93d3e9b7261/scratchpad"
IDS=[l.split("\t") for l in open(WORK+"/ids_1000.txt").read().splitlines() if l.strip()]
RESULTS=WORK+"/scale_results.jsonl"
BATCH=100

ALIAS={"MEK1":"MAP2K1","MEK2":"MAP2K2","ERK1":"MAPK3","ERK2":"MAPK1","4EBP1":"EIF4EBP1",
       "PUMA":"BBC3","S6K1":"RPS6KB1","TGFb1":"TGFB1","TGFb2":"TGFB2"}

import db_builder
from db_builder import BASELINE_NETWORK, DISEASE_SNT_SIGNATURES
panel=set()
for h,s,*_ in BASELINE_NETWORK: panel|={h,s}
for r in DISEASE_SNT_SIGNATURES: panel|={r[1],r[2]}
wanted={ALIAS.get(g,g):g for g in panel}

# build DB once
if not os.path.exists(DB):
    c=sqlite3.connect(DB); db_builder.create_schema(c); db_builder.seed_baseline(c)
    db_builder.seed_disease_signatures(c); db_builder.seed_healing_patterns(c); c.close()

import agent_logic as A

done=set()
if os.path.exists(RESULTS):
    for l in open(RESULTS):
        try: done.add(json.loads(l)["file_id"])
        except: pass

fid2case={fid:case for fid,case in IDS}
todo=[(fid,case) for fid,case in IDS if fid not in done]
print(f"[scale] total={len(IDS)} done={len(done)} todo={len(todo)}",flush=True)

def extract_panel(tsv_text):
    out={}
    for line in tsv_text.splitlines():
        if not line or line[0]=="#" or line.startswith(("gene_id","N_")): continue
        p=line.split("\t")
        if len(p)>=7 and p[1] in wanted:
            try: out[wanted[p[1]]]=float(p[6])
            except: pass
    return out

def run_patient(case, expr):
    conn=sqlite3.connect(DB); cur=conn.cursor()
    for g,v in expr.items():
        cur.execute("INSERT OR REPLACE INTO patient_expression (patient_id,gene_id,tpm_value) VALUES (?,?,?)",(case,g,v))
    conn.commit(); conn.close()
    e=A.load_patient_expression(case)
    triage=A.run_level1_triage(e)
    confirmed=[m for m in triage if m.confirmed]
    already={(m.hub_gene,m.satellite_gene) for m in triage}
    orphans=A.run_level2_block_scanner(e,already)
    collapsed=[m for m in triage if m.confirmed and m.expected_anomaly=="HUB_COLLAPSE"]
    aco=A.run_aco_analysis(case,collapsed,e)
    # clean up patient rows to keep DB small
    conn=sqlite3.connect(DB); conn.execute("DELETE FROM patient_expression WHERE patient_id=?",(case,)); conn.commit(); conn.close()
    return {"file_id":None,"case":case,"genes":len(e),"confirmed":len(confirmed),
            "orphans":len(orphans),"aco_hubs":len(aco),
            "diseases":sorted({m.disease_name for m in confirmed}),
            "modes":sorted({r.collapse_mode for r in aco})}

out=open(RESULTS,"a")
for bi in range(0,len(todo),BATCH):
    batch=todo[bi:bi+BATCH]
    ids=[fid for fid,_ in batch]
    t0=time.time()
    for attempt in range(4):
        try:
            req=urllib.request.Request("https://api.gdc.cancer.gov/data",
                data=json.dumps({"ids":ids}).encode(),headers={"Content-Type":"application/json"})
            blob=urllib.request.urlopen(req,timeout=300).read()
            break
        except Exception as e:
            if attempt==3: raise
            time.sleep(2**(attempt+1))
    tf=tarfile.open(fileobj=io.BytesIO(blob),mode="r:gz")
    n_ok=0
    for m in tf.getmembers():
        if not m.name.endswith(".tsv"): continue
        folder=m.name.split("/")[0]  # = file_id
        case=fid2case.get(folder,folder)
        try:
            expr=extract_panel(tf.extractfile(m).read().decode())
            rec=run_patient(case,expr); rec["file_id"]=folder
            out.write(json.dumps(rec)+"\n"); out.flush()
            n_ok+=1
        except Exception as e:
            out.write(json.dumps({"file_id":folder,"case":case,"error":str(e)})+"\n"); out.flush()
    print(f"[scale] batch {bi//BATCH+1}: {n_ok} patients in {time.time()-t0:.0f}s (cum done={len(done)+bi+n_ok})",flush=True)
    time.sleep(2)  # be gentle with the API
out.close()

# aggregate
recs=[json.loads(l) for l in open(RESULTS) if l.strip()]
ok=[r for r in recs if "error" not in r]
errs=[r for r in recs if "error" in r]
agg={"n_patients":len(ok),"n_errors":len(errs),
     "confirmed_mean":round(statistics.mean(r["confirmed"] for r in ok),2) if ok else 0,
     "confirmed_median":statistics.median(r["confirmed"] for r in ok) if ok else 0,
     "orphans_mean":round(statistics.mean(r["orphans"] for r in ok),2) if ok else 0,
     "aco_mean":round(statistics.mean(r["aco_hubs"] for r in ok),2) if ok else 0}
# disease frequency
from collections import Counter
dc=Counter()
for r in ok:
    for d in r["diseases"]: dc[d]+=1
agg["disease_freq"]={k:v for k,v in dc.most_common()}
json.dump(agg,open(WORK+"/scale_summary.json","w"),indent=2)
print("[scale] DONE",json.dumps(agg)[:300],flush=True)
