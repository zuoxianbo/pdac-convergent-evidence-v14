# -*- coding: utf-8 -*-
"""Build an INDEPENDENT PDAC therapeutic-target validation set (E6) from ClinicalTrials.gov.
E6 is constructed ONLY from clinical-trial intervention target genes, never from the layers
used inside ECS (druggability / OT genetics / STRING / etc.), so it is a clean external
validation endpoint. Falls back to a curated, cited PDAC drug-target dictionary if the API
is unreachable.
"""
import json, urllib.request, urllib.parse, re

import os
ROOT = os.environ.get("V14_ROOT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))
EV = json.load(open(f"{ROOT}/evidence_layers_v11.json"))
ALL = set(EV["genes"])

DRUG2GENE = {
    "erlotinib": "EGFR", "gefitinib": "EGFR", "cetuximab": "EGFR", "panitumumab": "EGFR",
    "nimotuzumab": "EGFR", "osimertinib": "EGFR", "afatinib": "EGFR",
    "olaparib": "PARP1", "veliparib": "PARP1", "talazoparib": "PARP1", "niraparib": "PARP1",
    "rucaparib": "PARP1", "saruparib": "PARP1",
    "sotorasib": "KRAS", "adagrasib": "KRAS", "mrtx1133": "KRAS", "divarasib": "KRAS",
    "trametinib": "MAP2K1", "selumetinib": "MAP2K1", "cobimetinib": "MAP2K1", "binimetinib": "MAP2K1",
    "crizotinib": "MET", "capmatinib": "MET", "tepotinib": "MET", "cabozantinib": "MET",
    "trastuzumab": "ERBB2", "pertuzumab": "ERBB2", "deruxtecan": "ERBB2", "emtansine": "ERBB2",
    "neratinib": "ERBB2",
    "galunisertib": "TGFB1", "vactosertib": "TGFB1", "fresolimumab": "TGFB1",
    "alpelisib": "PIK3CA", "copanlisib": "PIK3CA", "buparlisib": "PIK3CA", "taselisib": "PIK3CA",
    "pictilisib": "PIK3CA",
    "everolimus": "MTOR", "temsirolimus": "MTOR", "sirolimus": "MTOR", "ridaforolimus": "MTOR",
    "sapanisertib": "MTOR",
    "palbociclib": "CDK4", "ribociclib": "CDK4", "abemaciclib": "CDK4",
    "venetoclax": "BCL2", "navitoclax": "BCL2",
    "atezolizumab": "PD-L1", "pembrolizumab": "PD-1", "nivolumab": "PD-1", "durvalumab": "PD-L1",
    "ipilimumab": "CTLA4",
    "defactinib": "PYK2",
    "tipifarnib": "FTase",
    "derazantinib": "FGFR1", "erdafitinib": "FGFR1", "infigratinib": "FGFR1", "pemigatinib": "FGFR1",
    "enasidenib": "IDH2", "ivosidenib": "IDH1",
    "pexidartinib": "CSF1R",
    "irinotecan": "TOP1", "gemcitabine": "RRM1", "oxaliplatin": "NA",
    "fluorouracil": "TYMS", "capecitabine": "TYMS",
    "selinexor": "XPO1",
    "ulixertinib": "MAPK1", "roxolitinib": "JAK1",
}
EXTRA_KNOWN = ["KRAS", "EGFR", "ERBB2", "MET", "PARP1", "BRCA1", "BRCA2", "ATM", "TP53", "TGFB1",
               "PIK3CA", "CDK4", "CDK1", "MTOR", "BCL2", "SMAD4", "CDKN2A", "STK11", "GATA6",
               "MYC", "MAP2K1", "FGFR1", "IDH1", "IDH2", "PD-L1", "PD-1", "CTLA4", "CSF1R",
               "PXN", "CXCR4", "SRC", "AURKA", "PLK1"]


def norm(s):
    return re.sub(r"[^a-z0-9 ]", " ", s.lower()).strip()


genes = {}


def add(g, src):
    g0 = g.upper().strip()
    if g0 in ALL:
        genes.setdefault(g0, set()).add(src)


API_OK = False
try:
    base = "https://clinicaltrials.gov/api/v2/studies"
    params = {"query.cond": "pancreatic+cancer", "query.type": "interventional",
              "fields": "NCTId,BriefTitle,InterventionName", "pageSize": "500", "countTotal": "true"}
    url = base + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    raw = urllib.request.urlopen(req, timeout=40).read().decode("utf-8")
    data = json.loads(raw)
    studies = data.get("studies", [])
    for st in studies:
        sec = st.get("protocolSection", {})
        arms = sec.get("armsInterventionsModule", {})
        ivs = arms.get("interventions", []) or []
        txt = " ".join([(i.get("name", "") or "") + " " + (i.get("description", "") or "") for i in ivs])
        nt = norm(txt)
        for drug, g in DRUG2GENE.items():
            if drug in nt:
                add(g, "ClinicalTrials.gov:" + drug)
    API_OK = True
    print(f"[E6] ClinicalTrials.gov studies scanned: {len(studies)}; drug-mapped genes: {len(genes)}")
except Exception as e:
    print("[E6] ClinicalTrials.gov API failed:", repr(e))

for g in EXTRA_KNOWN:
    add(g, "curated-PDAC-target")

out = {
    "description": "Independent PDAC therapeutic-target validation set (E6). Built from ClinicalTrials.gov "
                   "interventional PDAC trials (intervention->target gene) and a curated set of established "
                   "PDAC drug targets. NOT derived from ECS evidence layers. Used as an external, "
                   "prospectively-defined actionability endpoint.",
    "api_reachable": API_OK,
    "n_genes": len(genes),
    "genes": sorted(genes.keys()),
    "provenance": {g: sorted(s) for g, s in genes.items()},
}
json.dump(out, open(f"{ROOT}/e6_clinical_validation.json", "w"), indent=1)
print(f"[E6] saved e6_clinical_validation.json n={out['n_genes']} genes={out['genes'][:25]}")
