"""GreenLight · FastAPI stub implementing the module contract (docs/MODULE_API.md).
Replace `fit` and `forecast` with the tumour–immune ODE (or any other model).
Run:  pip install fastapi uvicorn && uvicorn api.fastapi_stub:app --port 8788
Then register http://localhost:8788/module on the Integration page.
"""
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import Any, Dict
import math

app = FastAPI(title="GreenLight module stub")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])

META = {
    "meta": {"id": "tumor_immune_py", "nameKey": "Tumour–immune module (Python)", "owner": "Tumour–immune Aim", "version": "0.1", "kind": "remote"},
    "consumes": {"variables": ["car_copies", "tumor_burden"], "planKeys": ["infusion.car_dose_frac", "ld.flu_frac", "ld.cy_frac"], "signals": []},
    "provides": {
        "signals": ["car_expansion"], "variables": ["car_copies", "tumor_burden"],
        "params": [{"id": "rho", "labelKey": "param.rho", "log": True, "lo": 0.2, "hi": 3}, {"id": "T0", "labelKey": "param.T0", "log": True, "lo": 1, "hi": 1500}],
        "outcomes": [{"id": "p_response", "labelKey": "out.p_response", "kind": "prob", "direction": "higher"},
                     {"id": "peak_expansion", "labelKey": "out.peak_expansion", "kind": "value", "direction": "higher"}],
    },
}

class FitBody(BaseModel):
    snapshot: Dict[str, Any]
    ctx: Dict[str, Any]

class ForecastBody(BaseModel):
    state: Dict[str, Any]
    plan: Dict[str, Any]
    ctx: Dict[str, Any]

@app.get("/module/meta")
def meta():
    return META

@app.post("/module/fit")
def fit(body: FitBody):
    obs = body.snapshot.get("observations", [])
    # TODO: posterior over the tumour–immune parameters given obs up to ctx.asOfDay
    return {"moduleId": "tumor_immune_py", "asOfDay": body.ctx.get("asOfDay"), "nObs": len(obs), "ess": 100,
            "params": {"rho": {"median": 2.0, "lo": 1.5, "hi": 2.6}, "T0": {"median": 80, "lo": 50, "hi": 120}}}

@app.post("/module/forecast")
def forecast(body: ForecastBody):
    horizon = int(body.ctx.get("horizon", 120))
    grid = list(range(-6, horizon + 1))
    dose = body.plan["stages"]["infusion"].get("car_dose_frac", 1.0)
    # TODO: simulate the ODE under body.plan for posterior draws in body.state
    E = [0.0 if d < 0 else dose * 5.0 * math.exp(-((d - 9) / 6.0) ** 2) for d in grid]
    T = [80.0 * math.exp(-0.05 * max(d, 0)) for d in grid]
    q = lambda arr, f: [v * f for v in arr]
    return {"grid": grid,
            "variables": {"car_copies": {"q05": q(E, 0.4), "q25": q(E, 0.7), "q50": E, "q75": q(E, 1.4), "q95": q(E, 2.2)},
                          "tumor_burden": {"q05": q(T, 0.5), "q25": q(T, 0.8), "q50": T, "q75": q(T, 1.3), "q95": q(T, 1.8)}},
            "outcomes": {"p_response": {"value": 0.7}, "peak_expansion": {"value": max(E), "lo": 0.5 * max(E), "hi": 2 * max(E)}},
            "signals": {"car_expansion": {"grid": grid, "values": [e / 5.0 for e in E]}}}
