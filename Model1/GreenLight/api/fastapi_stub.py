"""GreenLight (Model1) · FastAPI stub of a remote module (docs/MODULE_API.md, docs/ML_API.md).
GET /card returns the model card, POST /predict takes a patient snapshot as of a day.
Replace `predict` with your model.
Run:  pip install fastapi uvicorn && uvicorn api.fastapi_stub:app --port 8788
Then add http://localhost:8788 on Researcher -> Interfaces.
"""
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from typing import Any, Dict
import math

app = FastAPI(title="GreenLight remote module stub")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])

CARD = {"id": "anc_trend_py", "version": "0.1", "type": "ml",
        "outputs": [{"id": "p_nonquick", "labelKey": "out.p_nonquick", "unit": "%", "direction": "lower"}]}

@app.get("/card")
def card() -> Dict[str, Any]:
    return CARD

@app.post("/predict")
def predict(snapshot: Dict[str, Any]) -> Dict[str, Any]:
    # toy rule: low last ANC and a falling trend raise the probability of a non-quick recovery
    anc = [o for o in snapshot.get("obs", {}).get("anc", []) if o[0] >= 0]
    if not anc:
        return {"p_nonquick": None}
    last = max(anc[-1][1], 0.01)
    trend = (math.log(max(anc[-1][1], 0.01)) - math.log(max(anc[0][1], 0.01))) / max(anc[-1][0] - anc[0][0], 1)
    eta = -0.5 - 1.2 * math.log(last) - 4.0 * trend
    return {"p_nonquick": 1 / (1 + math.exp(-eta)), "landmark": snapshot.get("asOf")}
