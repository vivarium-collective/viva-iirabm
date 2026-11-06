"""Reproduction drivers for Cockrell & An (2018).

Each ``*_data`` function runs the real IIRABM to produce one paper figure's
underlying data; each ``fig_*`` function renders an interactive Plotly figure
in the paper's visual language (red = death, blue = heal). Studies' ``run.py``
scripts call these so the heavy logic lives in one tested place.
"""
from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .interventions import FIG2_DIGITIZED, TARGETS, flatten
from .runtime import DEATH_THRESHOLD, run_batch, run_sim

# paper-faithful palette (also legible in dark mode)
C_DEATH = "#c0392b"   # red trajectories -> death
C_HEAL = "#2d6fb3"    # blue trajectories -> heal
C_ACCENT = "#5b2d8f"  # intervention bars (paper's navy/purple)
PLOT_BG = "rgba(0,0,0,0)"

_LAYOUT = dict(
    paper_bgcolor=PLOT_BG, plot_bgcolor=PLOT_BG,
    font=dict(family="Inter, system-ui, sans-serif", size=13),
    margin=dict(l=60, r=20, t=50, b=50),
)


def _style(fig, **kw):
    fig.update_layout(**_LAYOUT, **kw)
    fig.update_xaxes(showgrid=True, gridcolor="rgba(128,128,128,0.15)", zeroline=False)
    fig.update_yaxes(showgrid=True, gridcolor="rgba(128,128,128,0.15)", zeroline=False)
    return fig


# ----------------------------------------------------------------- Fig 6
def fig6_data(n_seeds=100, n_steps=12500, **sim_kw):
    runs = run_batch(range(1, n_seeds + 1), n_steps=n_steps, **sim_kw)
    deaths = sum(r["died"] for r in runs)
    return {"runs": runs, "n_seeds": n_seeds, "n_steps": n_steps,
            "deaths": deaths, "mortality": deaths / n_seeds}


def fig6_figure(data, stride=15):
    fig = go.Figure()
    for r in data["runs"]:
        oxy = r["oxygen_deficit"].astype(float).copy()
        # render the death threshold as the ceiling for dead patients
        oxy[oxy < 0] = DEATH_THRESHOLD
        x = np.arange(len(oxy))
        # downsample for a compact, visually-identical figure
        xs, ys = x[::stride], oxy[::stride]
        color = C_DEATH if r["died"] else C_HEAL
        fig.add_trace(go.Scatter(
            x=xs, y=ys, mode="lines",
            line=dict(color=color, width=0.8), opacity=0.5,
            hoverinfo="skip", showlegend=False))
    # legend proxies
    fig.add_trace(go.Scatter(x=[None], y=[None], mode="lines",
                  line=dict(color=C_DEATH, width=2), name="death"))
    fig.add_trace(go.Scatter(x=[None], y=[None], mode="lines",
                  line=dict(color=C_HEAL, width=2), name="heal"))
    mort = 100 * data["mortality"]
    _style(fig, title=f"IIRABM stochastic trajectories — {mort:.0f}% mortality "
                      f"({data['deaths']}/{data['n_seeds']} seeds)  ·  paper: ~82%",
           xaxis_title="Timestep (7 min)", yaxis_title="Oxygen deficit (system damage)",
           legend=dict(orientation="h", y=1.08, x=1, xanchor="right"))
    return fig


# ----------------------------------------------------------------- Fig 4
def fig4_data(n_seeds=60, n_steps=12500, intervention=FIG2_DIGITIZED, **sim_kw):
    """Baseline vs treated; split treated into responders (heal) / non-responders."""
    base = run_batch(range(1, n_seeds + 1), n_steps=n_steps, **sim_kw)
    treated = run_batch(range(1, n_seeds + 1), n_steps=n_steps,
                        intervention=intervention, **sim_kw)
    responders = [r for r in treated if not r["died"]]
    nonresp = [r for r in treated if r["died"]]
    return {"base": base, "treated": treated,
            "responders": responders, "nonresponders": nonresp, "n_steps": n_steps,
            "base_mortality": sum(r["died"] for r in base) / len(base),
            "treated_mortality": sum(r["died"] for r in treated) / len(treated)}


def _avg_series(runs, key):
    if not runs:
        return None
    m = min(len(r[key]) for r in runs)
    arr = np.array([np.where(r[key][:m] < 0, DEATH_THRESHOLD if key == "oxygen_deficit" else 0.0,
                             r[key][:m]) for r in runs], dtype=float)
    return arr.mean(axis=0)


def fig4_figure(data):
    fig = make_subplots(rows=1, cols=3, subplot_titles=(
        "Oxygen deficit", "Total GCSF", "Neutrophil population"))
    resp, nonr = data["responders"], data["nonresponders"]
    # pick a single representative non-responder (first death), responder average
    one_nr = nonr[0] if nonr else None
    panels = [
        ("oxygen_deficit", 1),
        (None, 2),   # GCSF handled below (cytokine)
        (None, 3),   # neutrophils (cell pop)
    ]
    def series(run, kind):
        if kind == "GCSF":
            s = run["cytokines"]["GCSF"]
        elif kind == "pmn":
            s = run["cell_populations"]["pmn"]
        else:
            s = run[kind]
        return np.where(s < 0, np.nan, s)

    def resp_avg(kind):
        if not resp:
            return None
        m = min(len(series(r, kind)) for r in resp)
        return np.nanmean(np.array([series(r, kind)[:m] for r in resp]), axis=0)

    for col, kind in [(1, "oxygen_deficit"), (2, "GCSF"), (3, "pmn")]:
        ra = resp_avg(kind)
        if ra is not None:
            fig.add_trace(go.Scatter(y=ra, mode="lines", line=dict(color=C_HEAL, width=2),
                          name="responder avg", showlegend=(col == 1)), row=1, col=col)
        if one_nr is not None:
            s = series(one_nr, kind)
            fig.add_trace(go.Scatter(y=s, mode="lines", line=dict(color=C_DEATH, width=2),
                          name="non-responder", showlegend=(col == 1)), row=1, col=col)
    _style(fig, title=f"Non-responder dynamics — treated mortality "
                      f"{100*data['treated_mortality']:.0f}% vs baseline "
                      f"{100*data['base_mortality']:.0f}%",
           legend=dict(orientation="h", y=1.12, x=1, xanchor="right"))
    for c in (1, 2, 3):
        fig.update_xaxes(title_text="Timestep", row=1, col=c)
    return fig


# ----------------------------------------------------------------- Fig 2 / Table 1
def table1_data(seeds_n=60, n_steps=12500, genome=FIG2_DIGITIZED, **sim_kw):
    """Mortality vs number of sequential interventions applied (Table 1 shape)."""
    seeds = range(1, seeds_n + 1)
    rows = []
    for k in [0, 1, 3, 5, len(genome)]:
        g = genome[:k] if k else None
        runs = run_batch(seeds, n_steps=n_steps, intervention=g, **sim_kw)
        mr = sum(r["died"] for r in runs) / len(runs)
        rows.append({"n_interventions": k, "mortality": mr})
    return {"rows": rows, "seeds_n": seeds_n}


def fig2_figure(genome=FIG2_DIGITIZED):
    """Fig-2-style: 8 stacked bar panels, one per stage, log2(multiplier) per target."""
    n = len(genome)
    fig = make_subplots(rows=n, cols=1, shared_xaxes=True,
                        subplot_titles=[f"{12*(i+1)} hrs post injury" for i in range(n)])
    for i, stage in enumerate(genome):
        log2v = [float(np.log2(v)) for v in stage]
        fig.add_trace(go.Bar(x=TARGETS, y=log2v, marker_color=C_ACCENT,
                      showlegend=False), row=i + 1, col=1)
        fig.update_yaxes(range=[-5, 5], row=i + 1, col=1)
    _style(fig, title="GA-identified intervention (digitized Fig 2) — "
                      "log₂ augmentation/inhibition per cytokine",
           height=140 * n)
    return fig


def table1_figure(data):
    rows = data["rows"]
    fig = go.Figure(go.Bar(
        x=[str(r["n_interventions"]) for r in rows],
        y=[100 * r["mortality"] for r in rows],
        marker_color=C_ACCENT,
        text=[f"{100*r['mortality']:.0f}%" for r in rows], textposition="outside"))
    _style(fig, title="Mortality vs number of sequential interventions",
           xaxis_title="# interventions applied", yaxis_title="Mortality (%)",
           yaxis_range=[0, 100])
    return fig


def ga_convergence_figure(history):
    fig = go.Figure(go.Scatter(y=history, x=list(range(len(history))),
                    mode="lines+markers", line=dict(color=C_ACCENT, width=2)))
    _style(fig, title="GA convergence — best fitness (oxygen deficit + infection)",
           xaxis_title="Generation", yaxis_title="Best fitness (lower = healthier)")
    return fig


# ----------------------------------------------------------------- Fig 3
def fig3_data(param_sets, seeds_n=20, n_steps=9000, genome=FIG2_DIGITIZED):
    """For each parameter set: base mortality and post-intervention mortality."""
    seeds = range(1, seeds_n + 1)
    out = []
    for ps in param_sets:
        base = run_batch(seeds, n_steps=n_steps, **ps)
        treated = run_batch(seeds, n_steps=n_steps, intervention=genome, **ps)
        out.append({
            "params": ps,
            "base_mr": sum(r["died"] for r in base) / seeds_n,
            "post_mr": sum(r["died"] for r in treated) / seeds_n,
        })
    return {"points": out, "seeds_n": seeds_n}


def fig3_figure(data):
    pts = data["points"]
    bx = [p["base_mr"] for p in pts]
    py = [p["post_mr"] for p in pts]
    improved = sum(1 for p in pts if p["post_mr"] < p["base_mr"])
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode="lines",
                  line=dict(color="rgba(128,128,128,0.6)", dash="dash"),
                  name="no change", hoverinfo="skip"))
    fig.add_trace(go.Scatter(
        x=bx, y=py, mode="markers",
        marker=dict(size=9, color=[C_HEAL if p["post_mr"] < p["base_mr"] else C_DEATH for p in pts],
                    line=dict(width=0.5, color="rgba(0,0,0,0.3)")),
        text=[str(p["params"]) for p in pts], name="parameter set"))
    _style(fig, title=f"Generalizability — {improved}/{len(pts)} parameter sets improved "
                      "(below diagonal = intervention beneficial)",
           xaxis_title="Base mortality rate", yaxis_title="Post-intervention mortality rate",
           xaxis_range=[0, 1], yaxis_range=[0, 1], legend=dict(y=0.02, x=0.98, xanchor="right"))
    return fig
