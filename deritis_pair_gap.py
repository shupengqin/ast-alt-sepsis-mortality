"""How far apart are the first AST and the first ALT in eICU?"""

import pandas as pd

from deritis_paper_analysis import connect

conn = connect("eicu")
tab = pd.read_sql_query(
    """
    WITH sepsis AS (
        SELECT DISTINCT patientunitstayid FROM diagnosis
        WHERE lower(diagnosisstring) LIKE '%sepsis%' OR lower(diagnosisstring) LIKE '%septic%'
    ),
    labs AS (
        SELECT l.patientunitstayid, lower(l.labname) AS labname, l.labresult, l.labresultoffset
        FROM lab l
        JOIN sepsis s ON s.patientunitstayid = l.patientunitstayid
        WHERE l.labresultoffset >= 0 AND l.labresultoffset < 1440
          AND l.labresult IS NOT NULL
          AND lower(l.labname) IN ('ast (sgot)', 'alt (sgpt)')
    ),
    firsts AS (
        SELECT DISTINCT ON (patientunitstayid, labname)
               patientunitstayid, labname, labresultoffset
        FROM labs
        ORDER BY patientunitstayid, labname, labresultoffset, labresult
    ),
    wide AS (
        SELECT patientunitstayid,
               MAX(labresultoffset) FILTER (WHERE labname = 'ast (sgot)') AS ast_off,
               MAX(labresultoffset) FILTER (WHERE labname = 'alt (sgpt)') AS alt_off
        FROM firsts
        GROUP BY patientunitstayid
    )
    SELECT ast_off, alt_off
    FROM wide
    WHERE ast_off IS NOT NULL AND alt_off IS NOT NULL
    """,
    conn,
)
conn.close()
gap = (tab["ast_off"] - tab["alt_off"]).abs()
n = len(gap)
print("n_both", n)
print("same_offset", int((gap == 0).sum()), f"{(gap == 0).mean():.3f}")
for minutes in (15, 60, 120, 360):
    print(f"within_{minutes}m", int((gap <= minutes).sum()), f"{(gap <= minutes).mean():.3f}")
print("median", float(gap.median()), "p90", float(gap.quantile(0.9)), "p99", float(gap.quantile(0.99)))
print("max", int(gap.max()))
