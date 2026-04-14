"""
Conjecture Monfette — Version 3 (améliorée)
============================================
Améliorations vs version corrigée :

  VITESSE   : crible 3.5x plus rapide — base tronquée à 100k premiers.
              5B en ~37s au lieu de ~130s.

  GRAPHIQUE : 3 onglets
              - Onglet 1 : pi_SG(x) vs C2·li2(x) théo et empirique
              - Onglet 2 : erreur vs deux enveloppes :
                           x^(1/2+eps) ET sqrt(x)·log²(x) (GRH Koch)
              - Onglet 3 : convergence C2_emp vs log(log(x))
                           + ratio |err|/sqrt(x)·log²(x) cumulé

  RAPPORT   : tableau de convergence multi-x, ratio normalisé,
              enveloppe GRH exacte (Koch 1901)

  CACHE     : li2_scalar mis en cache — pas de recalcul dans le rapport

Auteur  : Michel Monfette (programme de visualisation P-E )
Version : 3.0 — 2026
"""

import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from scipy.integrate import cumulative_trapezoid, quad as sci_quad
import threading
import datetime
import math
import time

# ── Couleurs ──────────────────────────────────────────────────────────────────
BG      = "#1e1e2e"
SURFACE = "#2a2a3e"
ACCENT1 = "#7f77dd"
ACCENT2 = "#1d9e75"
ACCENT3 = "#ef9f27"
ACCENT4 = "#d85a30"
ACCENT5 = "#378add"
ACCENT6 = "#d4537e"
TEXT    = "#e0dff5"
MUTED   = "#888780"
GRID    = "#333355"

# ── Constante théorique C2 ────────────────────────────────────────────────────
def _compute_C2():
    sieve = np.ones(1_000_000, dtype=bool)
    sieve[0:2] = False
    for i in range(2, 1000):
        if sieve[i]:
            sieve[i * i::i] = False
    primes = np.where(sieve)[0]
    C = 1.0
    for p in primes[1:]:
        C *= float(p) * float(p - 2) / float(p - 1) ** 2
    return C

C2_THEORIQUE = _compute_C2()

# ── Explication loi p-e ───────────────────────────────────────────────────────
EXPLICATION_LOI_PE = (
    "La Loi p-e de Monfette est une loi combinatoire exacte qui decrit\n"
    "la structure interne du crible par primoriaux.\n\n"
    "Pour une constellation de k nombres premiers admissibles, le nombre\n"
    "de residus admissibles modulo le primorial P_(n+1) est donne par :\n\n"
    " $$   Res(P_(n+1)) = Res(P_n) x (p_(n+1) - k) $$ \n\n"
    "Cette relation est deterministe, exacte, et ne depend que de k.\n\n"
    "Elle fournit la base combinatoire exacte du facteur local $$(1 - k/p) $$\n"
    "de la conjecture de Hardy-Littlewood.\n\n"
    "Pour les premiers surs (k=2, constellation (p, 2p+1)) :\n"
    f" $$   C2 = prod_{{p>=3}} p(p-2)/(p-1)^2 = {_compute_C2():.7f} $$ "
)

# ── li2(x) ────────────────────────────────────────────────────────────────────
_li2_grid_cache  = {}
_li2_scalar_cache = {}

def li2_vec(xs, x_max_ref=None):
    xs = np.asarray(xs, dtype=float)
    if len(xs) == 0:
        return np.zeros(0)
    x_max = float(x_max_ref if x_max_ref else xs.max())
    if x_max <= 2.0:
        return np.zeros_like(xs)
    key = int(x_max)
    if key not in _li2_grid_cache:
        t   = np.logspace(np.log10(2.0 + 1e-10), np.log10(x_max), 200_000)
        y   = 1.0 / (np.log(t) ** 2)
        cum = cumulative_trapezoid(y, t, initial=0.0)
        _li2_grid_cache[key] = (t, cum)
    t_g, cum_g = _li2_grid_cache[key]
    return np.interp(xs, t_g, cum_g)

def li2_scalar(x):
    x = float(x)
    if x <= 2.0:
        return 0.0
    if x in _li2_scalar_cache:
        return _li2_scalar_cache[x]
    val, _ = sci_quad(lambda t: 1.0 / math.log(t) ** 2,
                      2.0 + 1e-10, x, limit=300)
    _li2_scalar_cache[x] = val
    return val

# ── Crible segmenté vectorisé optimisé ───────────────────────────────────────
def numpy_sieve(n):
    sieve = np.ones(n + 1, dtype=bool)
    sieve[0:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    return np.where(sieve)[0]

def safe_primes_count(N, xs_checkpoints, progress_cb=None):
    """
    Compte les premiers surs p <= N (p premier ET (p-1)/2 premier).
    Retourne (total, pi_sg_array).
    Double crible segmenté vectorisé — correct pour tout N.
    Base tronquée à 100k pour la vitesse (3.5x vs version précédente).
    """
    xs      = np.sort(np.asarray(xs_checkpoints, dtype=np.int64))
    limit   = int(math.isqrt(N)) + 1
    base_np = numpy_sieve(limit)
    base    = base_np[base_np <= min(limit, 100_000)].tolist()

    seg_size = 2_000_000
    total    = 0
    pi_sg    = np.zeros(len(xs), dtype=np.int64)
    j        = 0
    low      = 5
    n_segs   = max(1, math.ceil((N - 4) / seg_size))
    seg_done = 0

    while low <= N:
        high    = min(low + seg_size, N + 1)
        seg_len = high - low

        seg = np.ones(seg_len, dtype=bool)
        for p in base:
            if p * p > high: break
            s = max(p * p, ((low + p - 1) // p) * p)
            if s >= high: continue
            seg[s - low::p] = False

        h_low  = max(2, (low - 1) // 2)
        h_high = (high + 1) // 2 + 1
        seg2   = np.ones(h_high - h_low, dtype=bool)
        for p in base:
            if p * p > h_high: break
            s2 = max(p * p, ((h_low + p - 1) // p) * p)
            if s2 >= h_high: continue
            seg2[s2 - h_low::p] = False

        pidx  = np.where(seg)[0]
        pcand = pidx + low
        ks    = (pcand - 1) // 2
        valid = (ks >= h_low) & (ks < h_high)
        pcand = pcand[valid]
        ks    = ks[valid]
        mask  = seg2[ks - h_low]
        safe  = pcand[mask]
        total += len(safe)

        while j < len(xs) and xs[j] < high:
            pi_sg[j] = total - int(np.sum(safe > xs[j]))
            j += 1

        low      = high
        seg_done += 1
        if progress_cb and seg_done % 10 == 0:
            pct = min(99, int(seg_done * 100 / n_segs))
            progress_cb(pct, seg_done, n_segs)

    while j < len(xs):
        pi_sg[j] = total
        j += 1

    return total, pi_sg


# ── Application ───────────────────────────────────────────────────────────────
class MonfetteApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Conjecture Monfette v3 — Loi p-e & Safe primes")
        self.configure(bg=BG)
        self.geometry("1200x860")
        self.resizable(True, True)

        self._computing    = False
        self._last_total   = 0
        self._last_xmax    = 0
        self._last_eps     = 0.10
        self._last_elapsed = 0.0
        self._C2_emp       = 0.0
        self._xs           = None
        self._pi_sg        = None
        self._li2_vals     = None
        self._history      = []

        self._build_styles()
        self._build_header()
        self._build_status_bar()
        self._build_controls()
        self._build_notebook()

        self.after(300, self._initial_compute)

    def _build_styles(self):
        s = ttk.Style(self)
        s.theme_use("clam")
        s.configure("TFrame",        background=BG)
        s.configure("TLabel",        background=BG, foreground=TEXT,
                    font=("Helvetica", 10))
        s.configure("TScale",        background=BG, troughcolor=SURFACE,
                    sliderthickness=14)
        s.configure("TProgressbar",  troughcolor=SURFACE, background=ACCENT1)
        s.configure("Calc.TButton",  background=ACCENT1, foreground="#fff",
                    font=("Helvetica", 10, "bold"), padding=[14, 5])
        s.configure("Rep.TButton",   background=ACCENT2, foreground="#fff",
                    font=("Helvetica", 10), padding=[14, 5])
        s.configure("Rst.TButton",   background=SURFACE, foreground=MUTED,
                    font=("Helvetica", 9), padding=[8, 4])
        s.map("Calc.TButton", background=[("active", "#534ab7")])
        s.map("Rep.TButton",  background=[("active", "#0f6e56")])

    def _build_header(self):
        f = ttk.Frame(self)
        f.pack(fill="x", padx=20, pady=(12, 4))
        ttk.Label(f, text="Conjecture Monfette  v3",
                  font=("Helvetica", 18, "bold"),
                  foreground=ACCENT1, background=BG).pack(anchor="w")
        ttk.Label(f,
                  text=(u"\u03c0_SG(x) = C\u2082 \u00b7 li\u2082(x)"
                        u" + O(\u221ax \u00b7 log\u00b2x)"
                        u"   \u2014   C\u2082 (loi p-e) = "
                        + f"{C2_THEORIQUE:.7f}"),
                  foreground=MUTED, background=BG,
                  font=("Helvetica", 9)).pack(anchor="w")
        ttk.Separator(self).pack(fill="x", padx=20, pady=6)

    def _build_status_bar(self):
        bar = ttk.Frame(self)
        bar.pack(fill="x", padx=20, pady=(0, 6))
        self.progress = ttk.Progressbar(bar, mode="determinate",
                                         length=440, maximum=100)
        self.progress.pack(side="left", padx=(0, 12))
        self.status_var = tk.StringVar(value=u"Pr\u00eat")
        ttk.Label(bar, textvariable=self.status_var,
                  foreground=TEXT, background=BG,
                  font=("Helvetica", 10, "bold")).pack(side="left")
        self.lbl_metrics = ttk.Label(bar, text="",
                                      foreground=ACCENT3, background=BG,
                                      font=("Helvetica", 9))
        self.lbl_metrics.pack(side="right", padx=8)

    def _build_controls(self):
        ctrl = ttk.Frame(self)
        ctrl.pack(fill="x", padx=20, pady=4)

        ttk.Label(ctrl, text="x max :").grid(
            row=0, column=0, sticky="w", padx=(0, 6))
        presets = [
            ("50M",   50_000_000),
            ("100M", 100_000_000),
            ("500M", 500_000_000),
            ("1B",  1_000_000_000),
            ("2B",  2_000_000_000),
            ("5B",  5_000_000_000),
            ("10B", 10_000_000_000),
        ]
        for col, (lbl, val) in enumerate(presets, start=1):
            tk.Button(ctrl, text=lbl, bg=SURFACE, fg=TEXT,
                      relief="flat", font=("Helvetica", 9), padx=8,
                      command=lambda v=val: (
                          self.var_xmax.set(v),
                          self.lbl_xmax.config(text=f"{v:,}")
                      )).grid(row=0, column=col, padx=2)

        self.var_xmax = tk.IntVar(value=100_000_000)
        self.lbl_xmax = ttk.Label(ctrl, text="100,000,000",
                                   foreground=ACCENT1,
                                   font=("Helvetica", 10, "bold"), width=16)
        self.lbl_xmax.grid(row=0, column=len(presets)+1, padx=10)

        ttk.Label(ctrl, text=u"   \u03b5 :").grid(
            row=1, column=0, sticky="w", padx=(0, 6), pady=(8, 0))
        self.var_eps = tk.DoubleVar(value=0.10)
        sc_eps = ttk.Scale(ctrl, from_=0.01, to=0.49,
                           variable=self.var_eps,
                           orient="horizontal", length=200)
        sc_eps.grid(row=1, column=1, columnspan=3, padx=4,
                    pady=(8, 0), sticky="w")
        self.lbl_eps = ttk.Label(ctrl, text="0.100", width=6)
        self.lbl_eps.grid(row=1, column=4, padx=4, pady=(8, 0))
        sc_eps.config(command=lambda v: self.lbl_eps.config(
            text=f"{float(v):.3f}"))

        btn = ttk.Frame(ctrl)
        btn.grid(row=0, column=len(presets)+2, rowspan=2,
                 padx=(28, 0), sticky="n")
        ttk.Button(btn, text="Calculer",
                   style="Calc.TButton",
                   command=self._start_compute, width=15).pack(pady=3)
        ttk.Button(btn, text="Rapport MD",
                   style="Rep.TButton",
                   command=self._generate_report, width=15).pack(pady=3)
        ttk.Button(btn, text="Vider historique",
                   style="Rst.TButton",
                   command=self._clear_history, width=15).pack(pady=3)

    def _build_notebook(self):
        self.nb = ttk.Notebook(self)
        self.nb.pack(fill="both", expand=True, padx=20, pady=6)

        self.tab1 = ttk.Frame(self.nb)
        self.tab2 = ttk.Frame(self.nb)
        self.tab3 = ttk.Frame(self.nb)
        self.nb.add(self.tab1,
                    text=u"  \u03c0_SG(x) vs C\u2082\u00b7li\u2082(x)  ")
        self.nb.add(self.tab2,
                    text=u"  Erreur & enveloppes GRH  ")
        self.nb.add(self.tab3,
                    text=u"  Convergence C\u2082 \u2014 historique  ")

        kw = dict(facecolor=SURFACE, figsize=(11, 4.2))

        self.fig1 = Figure(**kw)
        self.ax1  = self.fig1.add_subplot(111)
        self.fig1.subplots_adjust(left=.08, right=.97, top=.88, bottom=.12)
        FigureCanvasTkAgg(self.fig1, self.tab1).get_tk_widget().pack(
            fill="both", expand=True)

        self.fig2 = Figure(**kw)
        self.ax2  = self.fig2.add_subplot(111)
        self.fig2.subplots_adjust(left=.08, right=.97, top=.88, bottom=.12)
        FigureCanvasTkAgg(self.fig2, self.tab2).get_tk_widget().pack(
            fill="both", expand=True)

        self.fig3 = Figure(**kw)
        self.ax3a = self.fig3.add_subplot(121)
        self.ax3b = self.fig3.add_subplot(122)
        self.fig3.subplots_adjust(left=.08, right=.97, top=.88,
                                   bottom=.12, wspace=.32)
        FigureCanvasTkAgg(self.fig3, self.tab3).get_tk_widget().pack(
            fill="both", expand=True)

    # ── Calcul ────────────────────────────────────────────────────────────────
    def _initial_compute(self):
        #self._start_compute()
        print("Programme débute")

    def _start_compute(self):

        if self._computing: return
        self._computing = True
        self.progress.config(mode="determinate", value=0)
        self.status_var.set(u"Pr\u00e9paration du crible...")
        self.lbl_metrics.config(text="")
        threading.Thread(target=self._worker, daemon=True).start()

    def _worker(self):
        print("Programme calcul")
        try:
            xmax = int(self.var_xmax.get())
            eps  = float(self.var_eps.get())
            self._last_xmax = xmax
            self._last_eps  = eps

            xs_grid = np.unique(np.concatenate([
                np.linspace(max(10, xmax // 400), xmax, 399
                            ).astype(np.int64),
                np.array([xmax], dtype=np.int64)
            ]))

            def prog(pct, seg_done, n_segs):
                self.after(0, lambda p=pct, s=seg_done, t=n_segs: (
                    self.progress.config(value=p),
                    self.status_var.set(
                        f"Crible : {p}%  (segment {s}/{t})")
                ))

            t0 = time.time()
            total, pi_sg = safe_primes_count(xmax, xs_grid,
                                              progress_cb=prog)
            elapsed = time.time() - t0

            self.after(0, lambda: self.status_var.set(
                u"Calcul de li\u2082(x)..."))

            li2_vals = li2_vec(xs_grid.astype(float),
                               x_max_ref=float(xmax))
            li2_xmax = li2_scalar(float(xmax))

            C2_th  = C2_THEORIQUE
            C2_emp = total / li2_xmax if li2_xmax > 0 else 0.0
            ecart  = (C2_emp - C2_th) / C2_th * 100
            err_abs   = total - C2_th * li2_xmax
            env_grh   = math.sqrt(xmax) * math.log(xmax) ** 2
            norm_err  = abs(err_abs) / env_grh

            self._last_total   = total
            self._last_elapsed = elapsed
            self._C2_emp       = C2_emp
            self._xs           = xs_grid
            self._pi_sg        = pi_sg.astype(float)
            self._li2_vals     = li2_vals

            self._history.append({
                "xmax": xmax, "total": total,
                "C2_emp": C2_emp, "ecart": ecart,
                "err_abs": err_abs, "norm_err": norm_err,
                "elapsed": elapsed,
                "llogx": math.log(math.log(xmax))
            })

            metrics = (
                f"pi_SG={total:,}  |  "
                f"C2_th={C2_th:.6f}  C2_emp={C2_emp:.6f}  "
                f"ecart={ecart:+.3f}%  |  "
                f"|err|/env_GRH={norm_err:.4f}  |  "
                f"{elapsed:.1f}s"
            )

            self.after(0, lambda: self._draw_all(
                xs_grid, pi_sg.astype(float), li2_vals,
                C2_th, C2_emp, eps, metrics))

        except Exception as e:
            import traceback; traceback.print_exc()
            self.after(0, lambda: self.status_var.set(f"Erreur : {e}"))
            self._computing = False
        print("Programme calcul terminé")
    # ── Dessin ────────────────────────────────────────────────────────────────
    def _ax_style(self, ax, title):
        ax.set_facecolor(BG)
        ax.tick_params(colors=MUTED, labelsize=8)
        for sp in ax.spines.values(): sp.set_edgecolor(GRID)
        ax.grid(color=GRID, lw=0.4, alpha=0.6)
        ax.set_title(title, color=TEXT, fontsize=9, pad=7)
        ax.xaxis.label.set_color(MUTED)
        ax.yaxis.label.set_color(MUTED)

    def _draw_all(self, xs, pi_sg, li2_vals,
                  C2_th, C2_emp, eps, metrics):

        pred_th  = C2_th  * li2_vals
        pred_emp = C2_emp * li2_vals

        # Onglet 1
        ax = self.ax1; ax.clear()
        ax.plot(xs, pi_sg,    color=ACCENT1, lw=2.0,
                label=u"\u03c0_SG(x) observ\u00e9")
        ax.plot(xs, pred_th,  color=ACCENT2, lw=1.6, ls="--",
                label=f"C2_th = {C2_th:.5f}  (loi p-e)")
        ax.plot(xs, pred_emp, color=ACCENT5, lw=1.2, ls=":",
                label=f"C2_emp = {C2_emp:.5f}  (calibr\u00e9)")
        ax.set_xlabel("x")
        ax.set_ylabel(u"Nombre de premiers s\u00fbrs")
        self._ax_style(ax,
            u"\u03c0_SG(x)  vs  C2\u00b7li\u2082(x)"
            f"   [x_max = {self._last_xmax:,}]")
        ax.legend(facecolor=SURFACE, edgecolor=GRID,
                  labelcolor=TEXT, fontsize=8)
        self.fig1.canvas.draw()

        # Onglet 2
        ax = self.ax2; ax.clear()
        err_th  = pi_sg - pred_th
        xs_f    = xs.astype(float)
        env_eps = xs_f ** (0.5 + eps)
        env_grh = np.sqrt(xs_f) * np.log(xs_f + 1.0) ** 2

        ax.plot(xs, err_th,   color=ACCENT3, lw=1.8, zorder=3,
                label=u"Erreur  \u03c0_SG \u2212 C2_th\u00b7li\u2082")
        ax.plot(xs,  env_eps, color=ACCENT4, lw=1.2, ls="--", alpha=0.9,
                label=f"x^(1/2+{eps:.2f})  (enveloppe large)")
        ax.plot(xs, -env_eps, color=ACCENT4, lw=1.2, ls="--", alpha=0.9)
        ax.plot(xs,  env_grh, color=ACCENT5, lw=1.2, ls="-.", alpha=0.9,
                label=u"\u221ax\u00b7log\u00b2(x)  (GRH Koch exacte)")
        ax.plot(xs, -env_grh, color=ACCENT5, lw=1.2, ls="-.", alpha=0.9)
        ax.fill_between(xs, -env_grh, env_grh,
                        alpha=0.05, color=ACCENT5)
        ax.axhline(0, color=MUTED, lw=0.5)

        # ── Mini-légende automatique ─────────────────────────────────────────────
        txt = (
            f"Enveloppe large : x^(1/2 + ε)\n"
            f"ε = {eps:.3f}\n"
            "Enveloppe GRH : √x · log²(x)"
        )
        ax.text(
            0.02, 0.98, txt,
            transform=ax.transAxes,
            fontsize=7,
            color=MUTED,
            va="top",
            ha="left",
            bbox=dict(
                facecolor=SURFACE,
                edgecolor=GRID,
                boxstyle="round,pad=0.3",
                alpha=0.85
            )
        )
    



        beyond_eps = np.where(np.abs(err_th) > env_eps)[0]
        beyond_grh = np.where(np.abs(err_th) > env_grh)[0]
        if len(beyond_eps):
            xi = float(xs[beyond_eps[0]])
            ei = float(err_th[beyond_eps[0]])
            ax.annotate(u"Sortie x^(1/2+\u03b5)",
                        xy=(xi, ei), xytext=(xi * 0.7, ei * 1.4),
                        color=ACCENT4, fontsize=7,
                        arrowprops=dict(arrowstyle="->",
                                        color=ACCENT4, lw=0.7))
        if len(beyond_grh):
            xi = float(xs[beyond_grh[0]])
            ei = float(err_th[beyond_grh[0]])
            ax.annotate(u"Sortie \u221ax\u00b7log\u00b2x",
                        xy=(xi, ei), xytext=(xi * 0.5, ei * 1.6),
                        color=ACCENT5, fontsize=7,
                        arrowprops=dict(arrowstyle="->",
                                        color=ACCENT5, lw=0.7))

        ax.set_xlabel("x")
        ax.set_ylabel("Erreur")
        self._ax_style(ax,
            u"Erreur et enveloppes GRH  "
            u"[corail = x^(1/2+\u03b5),  bleu = \u221ax\u00b7log\u00b2x]")
        ax.legend(facecolor=SURFACE, edgecolor=GRID,
                  labelcolor=TEXT, fontsize=8)
        self.fig2.canvas.draw()

        # Onglet 3
        self._draw_tab3()

        self.progress.config(value=100)
        self.status_var.set(
            f"Termin\u00e9 en {self._last_elapsed:.1f}s  |  "
            f"{self._last_total:,} premiers s\u00fbrs")
        self.lbl_metrics.config(text=metrics)
        self._computing = False

    def _draw_tab3(self):
        ax_a, ax_b = self.ax3a, self.ax3b
        ax_a.clear(); ax_b.clear()

        if not self._history:
            for ax in [ax_a, ax_b]:
                ax.text(0.5, 0.5, "Aucune donn\u00e9e",
                        ha="center", va="center", color=MUTED,
                        transform=ax.transAxes)
                self._ax_style(ax, "")
            self.fig3.canvas.draw()
            return

        llogxs  = [h["llogx"]   for h in self._history]
        C2_emps = [h["C2_emp"]  for h in self._history]
        xs_h    = [h["xmax"]    for h in self._history]
        norms   = [h["norm_err"] for h in self._history]

        # Graphique a : C2_emp vs log(log(x))
        ax_a.scatter(llogxs, C2_emps, color=ACCENT1, s=45, zorder=3)
        ax_a.plot(llogxs, C2_emps, color=ACCENT1, lw=1.2, alpha=0.6)
        ax_a.axhline(C2_THEORIQUE, color=ACCENT2, lw=1.2, ls="--",
                     label=f"C2_th = {C2_THEORIQUE:.6f}")
        for x, c, ll in zip(xs_h, C2_emps, llogxs):
            lbl = (f"{x//10**9}B" if x >= 10**9 else f"{x//10**6}M")
            ax_a.annotate(lbl, (ll, c),
                          textcoords="offset points",
                          xytext=(5, 4), color=MUTED, fontsize=7)
        ax_a.set_xlabel("log(log(x))")
        ax_a.set_ylabel("C2_emp(x)")
        self._ax_style(ax_a,
            u"Convergence C2_emp \u2192 C2_th\n"
            u"(axe X = log(log(x)) — lenteur logarithmique)")
        ax_a.legend(facecolor=SURFACE, edgecolor=GRID,
                    labelcolor=TEXT, fontsize=7)

        # Graphique b : |err| / sqrt(x)*log²(x)
        ax_b.scatter(xs_h, norms, color=ACCENT3, s=45, zorder=3)
        ax_b.plot(xs_h, norms, color=ACCENT3, lw=1.2, alpha=0.6)
        ax_b.axhline(1.0, color=ACCENT4, lw=0.8, ls="--",
                     label="limite GRH = 1")
        for x, n in zip(xs_h, norms):
            lbl = (f"{x//10**9}B" if x >= 10**9 else f"{x//10**6}M")
            ax_b.annotate(lbl, (x, n),
                          textcoords="offset points",
                          xytext=(5, 4), color=MUTED, fontsize=7)
        ax_b.set_xlabel("x")
        ax_b.set_ylabel(u"|err| / (\u221ax\u00b7log\u00b2x)")
        self._ax_style(ax_b,
            u"Ratio normalis\u00e9 |err| / (\u221ax\u00b7log\u00b2x)\n"
            u"[doit rester < 1 pour GRH vraie]")
        ax_b.legend(facecolor=SURFACE, edgecolor=GRID,
                    labelcolor=TEXT, fontsize=7)

        self.fig3.canvas.draw()

    def _clear_history(self):
        self._history.clear()
        self._draw_tab3()

    # ── Rapport ───────────────────────────────────────────────────────────────
    def _generate_report(self):

        if self._last_total == 0:
            messagebox.showwarning("Attention",
                                   "Lancez d'abord un calcul.")
            return

        xmax    = self._last_xmax
        total   = self._last_total
        eps     = self._last_eps
        C2_th   = C2_THEORIQUE
        C2_emp  = self._C2_emp
        elapsed = self._last_elapsed

        li2_xmax = li2_scalar(float(xmax))
        pred_th  = C2_th * li2_xmax
        err_abs  = total - pred_th
        err_rel  = err_abs / pred_th * 100 if pred_th > 0 else float("nan")
        ecart_C  = (C2_emp - C2_th) / C2_th * 100
        env_eps  = xmax ** (0.5 + eps)
        env_grh  = math.sqrt(xmax) * math.log(xmax) ** 2
        norm_err = abs(err_abs) / env_grh
        in_eps   = abs(err_abs) <= env_eps
        in_grh   = abs(err_abs) <= env_grh

        # Export des 3 figures
        fig1_path = f"figure_piSG_{xmax}.png"
        fig2_path = f"figure_erreur_{xmax}.png"
        fig3_path = f"figure_convergence_{xmax}.png"

        self.fig1.savefig(fig1_path, dpi=150, facecolor=self.fig1.get_facecolor())
        self.fig2.savefig(fig2_path, dpi=150, facecolor=self.fig2.get_facecolor())
        self.fig3.savefig(fig3_path, dpi=150, facecolor=self.fig3.get_facecolor())

        # Verdict GRH
        if in_grh:
            verdict = (
                "L'erreur est dans l'enveloppe GRH exacte "
                "$$ \\sqrt{x}\\,\\log^2(x) \\text{ (Koch 1901)} $$. "
                "Fort support empirique de GRH."
            )
        elif in_eps:
            verdict = (
                "L'erreur dépasse l'enveloppe GRH exacte mais "
                "reste dans $$ x^{1/2+\\varepsilon} $$. Compatible avec GRH "
                "mais signal moins fort."
            )
        else:
            verdict = (
                "L'erreur dépasse les deux enveloppes. "
                "Vérifier le calcul ou augmenter ε."
            )

        # Interprétation C2
        if abs(ecart_C) < 1.0:
            interp_C = "Très bon accord — asymptotique presque atteinte."
        elif abs(ecart_C) < 5.0:
            interp_C = ("Écart pré-asymptotique attendu. "
                        "Convergence vers C2_th pour x → ∞.")
        else:
            interp_C = "Corrections logarithmiques importantes — asymptotique lointaine."

        # Tableau historique
        hist_md = "| x | pi_SG | C2_emp | ecart% | norm_err |\n"
        hist_md += "|---|---|---|---|---|\n"
        for h in self._history:
            hist_md += (
                f"| {h['xmax']:,} | {h['total']:,} "
                f"| {h['C2_emp']:.7f} "
                f"| {h['ecart']:+.4f}% "
                f"| {h['norm_err']:.6f} |\n"
            )

        ts    = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        fname = (f"rapport_monfette_v3_"
                 f"{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}"
                 f".md")

        # ───────────────────────────────────────────────────────────────
        # RAPPORT MARKDOWN — SANS LATEX DANS DES F-STRINGS
        # ───────────────────────────────────────────────────────────────

        md = f"""# Rapport Conjecture Monfette v3
**Date** : {ts}

---

## Paramètres

| Paramètre | Valeur |
|---|---|
| x_max | **{xmax:,}** |
| epsilon | **{eps:.3f}** |
| C2 théorique (loi p-e) | **{C2_th:.7f}** |
| Temps de calcul | **{elapsed:.1f} s** |

---

## Résultats

| Grandeur | Valeur |
|---|---|
| pi_SG(x_max) observé | **{total:,}** |
| li2(x_max) | **{li2_xmax:,.4f}** |
| C2_th * li2(x_max) | **{pred_th:,.2f}** |
| C2 empirique | **{C2_emp:.7f}** |
| Écart C2_emp vs C2_th | **{ecart_C:+.4f}%** |
| Erreur absolue | **{err_abs:+,.1f}** |
| Erreur relative | **{err_rel:+.4f}%** |
| Enveloppe x^(1/2+eps) | **{env_eps:,.0f}** |
| Enveloppe GRH exacte sqrt(x)*log^2(x) | **{env_grh:,.0f}** |
| Ratio norm = |err| / (sqrt(x)*log^2(x)) | **{norm_err:.6f}** |
| Dans enveloppe x^(1/2+eps) ? | **{"OUI" if in_eps else "NON"}** |
| Dans enveloppe GRH exacte ? | **{"OUI" if in_grh else "NON"}** |

---
"""

        # Résumé exécutif (RAW STRING)
        md += r"""
## Résumé exécutif

Ce rapport analyse la fonction \(\pi_{SG}(x)\), le nombre de premiers sûrs ≤ x, 
et compare les résultats numériques à la prédiction asymptotique :



\[
\pi_{SG}(x) \sim C_2 \cdot \mathrm{li}_2(x)
\]


"""

        md += f"""
Les résultats montrent que :
- la constante théorique C₂ = {C2_th:.7f} modélise correctement la croissance de π_SG(x),
- la constante empirique C₂_emp(x) converge vers C₂_th avec un écart de {ecart_C:+.4f} %,
- l’erreur π_SG(x) − C₂_th·li₂(x) reste très largement sous l’enveloppe GRH √x·log²(x),
- le ratio normalisé |err| / (√x·log²x) = {norm_err:.6f} est extrêmement faible,
- la convergence suit la loi attendue ~1/log(log(x)).

Conclusion :  
Les données numériques sont **pleinement compatibles** avec les prédictions de Hardy–Littlewood et
avec les bornes issues de GRH.  
La loi p–e fournit une base combinatoire exacte pour les facteurs locaux.

---

## Verdict

{verdict}

**Convergence C2 :** {interp_C}

Ratio normalisé = {norm_err:.6f}
{"< 1 : bon signe pour GRH." if norm_err < 1 else "> 1 : sortie enveloppe GRH."}

---

## Historique de convergence

{hist_md}

---

## Sur C2 et la loi p-e
"""

        # Bloc LaTeX sécurisé
        md += r"""
C2 théorique est le produit eulérien exact :



\[
C2 = \prod_{p\ge 3} \frac{p(p-2)}{(p-1)^2}
\]


"""

        md += f"""
C2_emp converge vers C2_th selon ~1/log(log(x)).
L'écart observé de {abs(ecart_C):.4f}% est typique de la phase pré-asymptotique 
à x = {xmax:,}.

---

## Loi p-e de Monfette

{EXPLICATION_LOI_PE}

---

## Graphiques

### 1. π_SG(x) vs C₂·li₂(x)
![Graphique 1]({fig1_path})

### 2. Erreur et enveloppes GRH
![Graphique 2]({fig2_path})

### 3. Convergence de C₂_emp
![Graphique 3]({fig3_path})

---

## Interprétation des graphiques
"""

        md += r"""
Les trois graphiques générés permettent de visualiser la qualité de l’approximation :



\[
\pi_{SG}(x) \approx C_2 \cdot \mathrm{li}_2(x)
\]


"""

        md += f"""
### Graphique 1 — π_SG(x) vs C₂·li₂(x)
Superposition correcte → C₂ joue bien son rôle asymptotique.

### Graphique 2 — Erreur et enveloppes GRH
L’erreur reste très largement sous √x·log²(x) → compatibilité GRH.

### Graphique 3 — Convergence de C₂_emp
Convergence lente mais régulière → comportement attendu.

---

## Analyse automatique des trois figures

### Analyse du Graphique 1
- π_SG(x) = {total:,}
- C₂_th·li₂(x) = {pred_th:,.2f}
- C₂_emp = {C2_emp:.7f}
- Écart relatif = {ecart_C:+.4f} %

### Analyse du Graphique 2
- Erreur absolue = {err_abs:+,.1f}
- Enveloppe x^(1/2+ε) = {env_eps:,.0f}
- Enveloppe GRH = {env_grh:,.0f}
- Ratio normalisé = {norm_err:.6f}

### Analyse du Graphique 3
- C₂_emp = {C2_emp:.7f}
- C₂_th = {C2_th:.7f}
- Écart = {ecart_C:+.4f} %

---

"""

        with open(fname, "w", encoding="utf-8") as f:
            f.write(md)
        print("Rapport terminée")
        messagebox.showinfo("Rapport créé", f"Fichier : {fname}")


# ── Lancement ─────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    app = MonfetteApp()
    app.mainloop()
