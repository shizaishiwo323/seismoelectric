from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from manimlib import *


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "hyperframes" / "schakel2011_code_video" / "assets"
TIMELINE_PATH = ASSETS / "timeline.json"
DATA_PATH = ROOT / "tmp" / "schakel2011_code_video" / "pe_result_summary.json"

FONT = "Microsoft YaHei"
INK = "#263238"
MUTED = "#607D8B"
PAPER = "#F7F9FC"
BLUE = "#2E86AB"
TEAL = "#0EAD69"
ORANGE = "#F28C28"
RED = "#C44536"
PURPLE = "#7B61FF"
YELLOW = "#F7C948"
GRAY = "#9AA6B2"


def fit(mob, width=12.0, height=6.2):
    if mob.get_width() > width:
        mob.set_width(width)
    if mob.get_height() > height:
        mob.set_height(height)
    return mob


def txt(s, size=28, color=INK, width=12.0):
    mob = Text(s, font=FONT, font_size=size)
    mob.set_color(color)
    return fit(mob, width=width)


def mtx(s, size=38, color=INK, width=12.0):
    mob = Tex(s, font_size=size)
    mob.set_color(color)
    return fit(mob, width=width)


def box(label, color=BLUE, width=2.4, height=0.7, size=22):
    rect = RoundedRectangle(width=width, height=height, corner_radius=0.08)
    rect.set_fill(color, opacity=0.14).set_stroke(color, width=2)
    label_mob = txt(label, size=size, width=width - 0.18)
    return VGroup(rect, label_mob)


def math_box(formula, color=BLUE, width=1.4, height=0.7, size=30):
    rect = RoundedRectangle(width=width, height=height, corner_radius=0.08)
    rect.set_fill(color, opacity=0.14).set_stroke(color, width=2)
    label_mob = mtx(formula, size=size, width=width - 0.18)
    return VGroup(rect, label_mob)


def formula_card(formula, source, color=BLUE, width=5.4, height=1.35, fsize=32, csize=18):
    card = RoundedRectangle(width=width, height=height, corner_radius=0.1)
    card.set_fill(WHITE, opacity=0.96).set_stroke(color, width=2)
    f = mtx(formula, size=fsize, width=width - 0.32)
    c = txt(source, size=csize, color=MUTED, width=width - 0.32)
    group = VGroup(f, c).arrange(DOWN, buff=0.10)
    return VGroup(card, group)


def arrow_between(a, b, color=MUTED):
    return Arrow(a.get_right() + 0.06 * RIGHT, b.get_left() + 0.06 * LEFT, buff=0.05, color=color)


class Schakel2011CodeStory(Scene):
    CONFIG = {
        "camera_config": {"background_color": PAPER, "pixel_width": 1280, "pixel_height": 720},
    }

    def setup(self):
        if TIMELINE_PATH.exists():
            data = json.loads(TIMELINE_PATH.read_text(encoding="utf-8"))
            self.timeline = {item["key"]: item for item in data["segments"]}
            self.total_duration = float(data["total_duration"])
        else:
            keys = [
                "title", "story", "source_map", "time_phase", "rt_inputs", "dynamic",
                "epsilon_sigma", "modes", "interface", "normalization_source",
                "sommerfeld", "porous_policy", "time_domain", "pe_results",
                "waveform", "liu_context", "appendix_material",
                "appendix_modes_boundary", "appendix_source_diagnostics", "ending",
            ]
            self.timeline = {key: {"end": (i + 1) * 24.0} for i, key in enumerate(keys)}
            self.total_duration = len(keys) * 24.0
        self.active = VGroup()
        self.pe_data = json.loads(DATA_PATH.read_text(encoding="utf-8")) if DATA_PATH.exists() else {}

    def hold(self, key):
        target = float(self.timeline[key]["end"])
        self.wait(max(0.2, target - self.time))

    def clear_active(self):
        if len(self.active) > 0:
            self.play(FadeOut(self.active), run_time=0.5)
        self.active = VGroup()

    def scene_title(self, title, subtitle=None):
        t = txt(title, size=34, width=11.5).to_edge(UP, buff=0.28)
        line = Line(LEFT * 5.7, RIGHT * 5.7, color=BLUE).next_to(t, DOWN, buff=0.10)
        group = VGroup(t, line)
        self.play(FadeIn(t, shift=0.12 * DOWN), ShowCreation(line), run_time=0.7)
        if subtitle:
            sub = txt(subtitle, size=20, color=MUTED, width=11.3).next_to(line, DOWN, buff=0.10)
            self.play(FadeIn(sub), run_time=0.35)
            group.add(sub)
        return group

    def pore_network(self, width=4.5):
        matrix = RoundedRectangle(width=width, height=2.85, corner_radius=0.08)
        matrix.set_fill("#E6EEF5", opacity=1).set_stroke("#9AB4C8", width=2)
        centers = [(-1.7, 0.65), (-0.55, 0.2), (0.55, 0.72), (1.55, -0.2), (-1.25, -0.75), (0.1, -0.78)]
        pores = VGroup()
        for x, y in centers:
            pores.add(Circle(radius=0.34).set_fill(WHITE, 1).set_stroke(BLUE, 2).move_to(matrix.get_center() + x * RIGHT + y * UP))
        channels = VGroup()
        for a, b in [(0, 1), (1, 2), (2, 3), (1, 5), (4, 5)]:
            channels.add(Line(pores[a].get_center(), pores[b].get_center(), color=BLUE, stroke_width=9).set_opacity(0.32))
        charges = VGroup(
            txt("+", 24, RED).move_to(pores[0].get_center()),
            txt("-", 24, TEAL).move_to(pores[1].get_center()),
            txt("+", 24, RED).move_to(pores[2].get_center()),
            txt("-", 24, TEAL).move_to(pores[3].get_center()),
        )
        return VGroup(matrix, channels, pores, charges)

    def wave_panel(self, width=4.8, height=2.0):
        axes = Axes(
            x_range=(0, 6, 1),
            y_range=(-1.15, 1.15, 1),
            width=width,
            height=height,
            axis_config={"stroke_color": MUTED, "stroke_width": 2},
        )
        c1 = axes.get_graph(lambda x: 0.78 * np.exp(-0.42 * (x - 3.2) ** 2) * np.sin(5.2 * x), color=BLUE)
        c2 = axes.get_graph(lambda x: -0.55 * np.exp(-0.50 * (x - 3.22) ** 2) * np.sin(5.2 * x), color=RED)
        t0 = DashedLine(axes.c2p(2.7, -1.0), axes.c2p(2.7, 1.0), color=ORANGE)
        label = mtx(r"T_0", size=28, color=ORANGE).next_to(t0, UP, buff=0.03)
        return VGroup(axes, c1, c2, t0, label)

    def data_curve(self, axes, samples, ycol, color, ylog=True):
        pts = []
        for row in samples:
            if ycol not in row or row[ycol] in ("", None):
                continue
            x = float(row["Time_min"])
            y = float(row[ycol])
            if not np.isfinite(x) or not np.isfinite(y) or y <= 0:
                continue
            yy = math.log10(y) if ylog else y
            pts.append(axes.c2p(x, yy))
        if len(pts) < 2:
            return VGroup()
        curve = VMobject()
        curve.set_points_as_corners(pts)
        curve.set_stroke(color, width=4)
        dots = VGroup(*[Dot(p, radius=0.035, color=color) for p in pts[:: max(1, len(pts) // 6)]])
        return VGroup(curve, dots)

    def construct(self):
        bg = Rectangle(width=FRAME_WIDTH + 1, height=FRAME_HEIGHT + 1)
        bg.set_fill(PAPER, opacity=1).set_stroke(width=0)
        self.add(bg)
        self.title_segment()
        self.story_segment()
        self.source_map_segment()
        self.time_phase_segment()
        self.rt_inputs_segment()
        self.dynamic_segment()
        self.epsilon_sigma_segment()
        self.modes_segment()
        self.interface_segment()
        self.normalization_source_segment()
        self.sommerfeld_segment()
        self.porous_policy_segment()
        self.time_domain_segment()
        self.pe_results_segment()
        self.waveform_segment()
        self.liu_context_segment()
        self.appendix_material_segment()
        self.appendix_modes_boundary_segment()
        self.appendix_source_diagnostics_segment()
        self.ending_segment()
        self.wait(max(0.2, self.total_duration - self.time))

    def title_segment(self):
        self.clear_active()
        title = txt("Schakel 2011 零偏移震电正演：从代码公式到研究故事", 39, width=11.4)
        chain = mtx(r"\{\phi,k_0,\alpha_\infty,c_H\}\rightarrow\{L,\sigma,\bar{\varepsilon}\}\rightarrow R^E,T^{TM}\rightarrow u(z,t)", 36, BLUE, width=11.0)
        claim = txt("溶蚀增强水力连通性，但 interface EM 峰值跨数量级降低", 25, RED, width=10.5)
        group = VGroup(title, chain, claim).arrange(DOWN, buff=0.35).shift(UP * 0.65)
        pore = self.pore_network(4.1).scale(0.72).to_edge(LEFT, buff=0.7).shift(DOWN * 1.25)
        wave = self.wave_panel(4.4, 1.85).to_edge(RIGHT, buff=0.75).shift(DOWN * 1.25)
        arrow = Arrow(pore.get_right(), wave.get_left(), color=ORANGE, stroke_width=6)
        self.play(FadeIn(group), run_time=1.0)
        self.play(FadeIn(pore), GrowArrow(arrow), FadeIn(wave[0]), ShowCreation(wave[1]), ShowCreation(wave[2]), FadeIn(wave[3:]), run_time=1.5)
        self.active = VGroup(group, pore, arrow, wave)
        self.hold("title")

    def story_segment(self):
        self.clear_active()
        header = self.scene_title("先确定故事主角", "公式不是清单，而是解释 interface EM 为什么变弱的证据链")
        main = box("主角：interface EM response", RED, width=4.4, height=0.95, size=26).move_to(UP * 1.1)
        actors = VGroup(
            box("孔隙率/渗透率", TEAL, width=2.5),
            box("H+ / zeta", TEAL, width=2.2),
            box("动态耦合/电导", BLUE, width=2.6),
            box("界面六阶系统", ORANGE, width=2.6),
            box("Sommerfeld 积分", PURPLE, width=2.8),
        ).arrange(RIGHT, buff=0.16).next_to(main, DOWN, buff=0.75)
        arrows = VGroup(*[Arrow(actors[i].get_top(), main.get_bottom(), color=MUTED, buff=0.08) for i in range(len(actors))])
        takeaway = txt("研究结论线：hydraulic connectivity ↑  不等于  electrokinetic imbalance ↑", 26, RED, width=10.0).to_edge(DOWN, buff=0.65)
        self.play(FadeIn(header), FadeIn(main), run_time=0.9)
        self.play(FadeIn(actors, lag_ratio=0.10), ShowCreation(arrows, lag_ratio=0.08), run_time=1.4)
        self.play(FadeIn(takeaway), run_time=0.6)
        self.active = VGroup(header, main, actors, arrows, takeaway)
        self.hold("story")

    def source_map_segment(self):
        self.clear_active()
        header = self.scene_title("代码公式的论文来源地图", "每个公式卡都要知道自己从哪篇论文来")
        left = VGroup(
            box("Schakel & Smeulders 2010 JASA", BLUE, width=4.2),
            box("动态系数 + 孔弹波模 + 界面矩阵", BLUE, width=4.6, size=20),
        ).arrange(DOWN, buff=0.18).to_edge(LEFT, buff=0.65).shift(UP * 0.8)
        mid = VGroup(
            box("Schakel et al. 2011 JAP/Geophysics", ORANGE, width=4.6),
            box("压力源 + Sommerfeld + 实验解释", ORANGE, width=4.6, size=20),
        ).arrange(DOWN, buff=0.18).move_to(UP * 0.8)
        right = VGroup(
            box("Liu et al. 2018 finite-offset VSEP", PURPLE, width=4.5),
            box("有限偏移频率-波数背景", PURPLE, width=4.1, size=20),
        ).arrange(DOWN, buff=0.18).to_edge(RIGHT, buff=0.65).shift(UP * 0.8)
        code = box("schakel2011_sommerfeld.py：offset_D = 0，主核为 Schakel 2011 Sommerfeld", RED, width=9.6, height=0.9, size=23).to_edge(DOWN, buff=0.75)
        arrows = VGroup(Arrow(left.get_bottom(), code.get_top(), color=BLUE), Arrow(mid.get_bottom(), code.get_top(), color=ORANGE), Arrow(right.get_bottom(), code.get_top(), color=PURPLE))
        self.play(FadeIn(header), FadeIn(left), FadeIn(mid), FadeIn(right), run_time=1.2)
        self.play(ShowCreation(arrows, lag_ratio=0.10), FadeIn(code), run_time=1.2)
        self.active = VGroup(header, left, mid, right, arrows, code)
        self.hold("source_map")

    def time_phase_segment(self):
        self.clear_active()
        header = self.scene_title("先把两条时间轴和相位说清", "这是后面判断 T0 和因果性的地基")
        slow = NumberLine(x_range=(0, 10, 2), width=5.0, color=TEAL).shift(UP * 1.2 + LEFT * 2.4)
        fast = NumberLine(x_range=(0, 10, 2), width=5.0, color=BLUE).shift(DOWN * 0.95 + RIGHT * 2.35)
        labels = VGroup(
            txt("dissolution time / Time_s：秒到小时", 25, TEAL, width=5.5).next_to(slow, UP),
            txt("waveform time t：微秒级传播", 25, BLUE, width=5.2).next_to(fast, DOWN),
            mtx(r"e^{i\omega t},\quad e^{-i\mathbf{k}\cdot\mathbf{x}}", 34, ORANGE, width=5.2).move_to(DOWN * 2.55 + LEFT * 2.4),
            mtx(r"z<0:\mathrm{fluid},\quad z>0:\mathrm{porous}", 32, PURPLE, width=5.2).move_to(DOWN * 2.55 + RIGHT * 2.8),
        )
        connect = Arrow(slow.get_right(), fast.get_left(), color=ORANGE)
        self.play(FadeIn(header), ShowCreation(slow), ShowCreation(fast), FadeIn(labels[:2]))
        self.play(GrowArrow(connect), FadeIn(labels[2:]), run_time=1.1)
        self.active = VGroup(header, slow, fast, labels, connect)
        self.hold("time_phase")

    def rt_inputs_segment(self):
        self.clear_active()
        header = self.scene_title("反应输运变量如何进入 Schakel 模型", "这一步是 RT-SE 映射，不是论文原始实验输入")
        inputs = VGroup(
            math_box(r"\phi", TEAL, 1.0),
            math_box(r"k_0", TEAL, 1.0),
            math_box(r"\alpha_\infty", TEAL, 1.15),
            math_box(r"c_H", TEAL, 1.0),
        ).arrange(DOWN, buff=0.15).to_edge(LEFT, buff=0.85)
        formulas = VGroup(
            formula_card(r"k_0[\mathrm{m^2}]=k_0[\mathrm{mD}]\times9.869233\times10^{-16}", "code unit conversion", BLUE, width=6.5),
            formula_card(r"pH=-\log_{10}\max(c_H,10^{-7})", "code RT interface", BLUE, width=6.5),
            formula_card(r"\zeta=\left[0.010+0.025\log_{10}C\right]{pH-2\over5}", "Schakel & Smeulders 2010 Appendix A, Eq. A5", ORANGE, width=6.5),
        ).arrange(DOWN, buff=0.22).to_edge(RIGHT, buff=0.6).shift(DOWN * 0.08)
        arrows = VGroup(*[Arrow(inputs.get_right(), formulas[i].get_left(), color=MUTED, buff=0.15) for i in range(3)])
        self.play(FadeIn(header), FadeIn(inputs, lag_ratio=0.12), run_time=0.8)
        self.play(ShowCreation(arrows, lag_ratio=0.12), FadeIn(formulas, lag_ratio=0.12), run_time=1.5)
        self.active = VGroup(header, inputs, formulas, arrows)
        self.hold("rt_inputs")

    def dynamic_segment(self):
        self.clear_active()
        header = self.scene_title("动态渗透率、Debye 长度、耦合系数和电导率", "Schakel & Smeulders 2010 Appendix A 是代码的材料参数发动机")
        cards = VGroup(
            formula_card(r"k(\omega)=k_0\left[\sqrt{1+i{\omega M\over2\omega_t}}+i{\omega\over\omega_t}\right]^{-1}", "Eq. A1", BLUE, width=6.1),
            formula_card(r"\omega_t={\phi\eta\over\alpha_\infty k_0\rho_f},\quad \Lambda=\sqrt{8\alpha_\infty k_0/(\phi M)}", "Eq. A2-A3", TEAL, width=6.1),
            formula_card(r"d=\left[\sum_l{(ez_l)^2N_l\over\varepsilon_0\varepsilon_f k_BT}\right]^{-1/2}", "Eq. A11", PURPLE, width=6.1),
            formula_card(r"L(\omega)\propto-{\phi\over\alpha_\infty}{\varepsilon_0\varepsilon_f\zeta\over\eta}", "Eq. A4, simplified front factor", ORANGE, width=6.1),
        ).arrange(DOWN, buff=0.16).to_edge(LEFT, buff=0.55).shift(DOWN * 0.1)
        pore = self.pore_network(3.6).scale(0.9).to_edge(RIGHT, buff=0.75).shift(UP * 0.2)
        sigma = formula_card(r"\sigma(\omega)={\phi\sigma_f\over\alpha_\infty}\left[1+{2(C_{em}+C_{os})\over\sigma_f\Lambda}\right]", "Eq. A6-A10", RED, width=4.9).next_to(pore, DOWN, buff=0.35)
        ring = Circle(radius=0.78, color=ORANGE).move_to(pore[2][1].get_center()).set_stroke(ORANGE, width=4)
        self.play(FadeIn(header), FadeIn(cards, lag_ratio=0.08), run_time=1.5)
        self.play(FadeIn(pore), ShowCreation(ring), FadeIn(sigma), run_time=1.1)
        self.active = VGroup(header, cards, pore, ring, sigma)
        self.hold("dynamic")

    def epsilon_sigma_segment(self):
        self.clear_active()
        header = self.scene_title("不要混淆两个电导率", "多孔体动态电导率和上覆自由流体电导率进入不同方程")
        left = formula_card(r"\bar{\varepsilon}=\varepsilon-i{\sigma(\omega)\over\omega}+i{\eta L^2\over\omega k(\omega)}", "Schakel & Smeulders 2010 Eq. 20", BLUE, width=5.7).to_edge(LEFT, buff=0.75).shift(UP * 0.8)
        right = formula_card(r"s_E^2=\mu_0\varepsilon_0\varepsilon_{fl}-i{\mu_0\sigma_{fl}\over\omega}", "fluid-side EM slowness", ORANGE, width=5.7).to_edge(RIGHT, buff=0.75).shift(UP * 0.8)
        note = VGroup(
            box("多孔介质动态体电导率", BLUE, width=4.1),
            box("上覆自由流体电导率", ORANGE, width=4.1),
            box("默认 constant，dynamic_pore_fluid 是敏感性假设", RED, width=6.2),
        ).arrange(DOWN, buff=0.22).move_to(DOWN * 1.25)
        neq = mtx(r"\sigma(\omega)\ne\sigma_{fl}", 42, RED).move_to(ORIGIN + UP * 0.2)
        self.play(FadeIn(header), FadeIn(left), FadeIn(right), run_time=1.0)
        self.play(FadeIn(neq), FadeIn(note, lag_ratio=0.10), run_time=1.1)
        self.active = VGroup(header, left, right, neq, note)
        self.hold("epsilon_sigma")

    def modes_segment(self):
        self.clear_active()
        header = self.scene_title("孔弹波模：把机械势和电势绑定", "从 Biot 系数到流固幅值比和电场幅值比")
        left = VGroup(
            formula_card(r"A,Q,R,\quad P=A+2G", "Schakel & Smeulders 2010 Eq. 8-10", TEAL, width=4.8),
            formula_card(r"\bar{\rho}_{11}=\rho_{11}-E_K,\ \bar{\rho}_{12}=\rho_{12}+E_K,\ \bar{\rho}_{22}=\rho_{22}-E_K", "Eq. 26-29", ORANGE, width=4.8),
        ).arrange(DOWN, buff=0.25).to_edge(LEFT, buff=0.55)
        modes = VGroup(box("fast P", BLUE, 1.45), box("slow P", BLUE, 1.45), box("TM", PURPLE, 1.25), box("SV", PURPLE, 1.25)).arrange_in_grid(2, 2, buff=0.25).move_to(ORIGIN)
        right = VGroup(
            formula_card(r"s_l^2={-d_1\over2d_2}\mp{1\over2}\sqrt{\left({d_1\over d_2}\right)^2-4{d_0\over d_2}}", "Eq. 24-25, 30", BLUE, width=4.9),
            formula_card(r"\beta={\mathrm{fluid\ potential}\over\mathrm{solid\ potential}},\quad \alpha={\mathrm{electric\ potential}\over\mathrm{mechanical\ potential}}", "Eq. 36-39", PURPLE, width=4.9),
        ).arrange(DOWN, buff=0.25).to_edge(RIGHT, buff=0.55)
        arrows = VGroup(Arrow(left.get_right(), modes.get_left(), color=MUTED), Arrow(modes.get_right(), right.get_left(), color=MUTED))
        self.play(FadeIn(header), FadeIn(left, lag_ratio=0.10), run_time=1.1)
        self.play(GrowArrow(arrows[0]), FadeIn(modes, lag_ratio=0.12), run_time=0.9)
        self.play(GrowArrow(arrows[1]), FadeIn(right, lag_ratio=0.10), run_time=1.1)
        self.active = VGroup(header, left, modes, right, arrows)
        self.hold("modes")

    def interface_segment(self):
        self.clear_active()
        header = self.scene_title("界面边界条件：六个未知量，六个方程", "Schakel & Smeulders 2010 Section III + Appendix B")
        water = Rectangle(width=10.8, height=1.75).set_fill("#DDEFF8", 1).set_stroke(width=0).to_edge(UP, buff=1.35)
        porous = Rectangle(width=10.8, height=1.75).set_fill("#EDE6D6", 1).set_stroke(width=0).to_edge(DOWN, buff=0.45)
        interface = Line(LEFT * 5.4, RIGHT * 5.4, color=INK, stroke_width=4).move_to(ORIGIN + DOWN * 0.1)
        inc = Arrow(LEFT * 4.0 + UP * 1.25, LEFT * 1.2 + DOWN * 0.08, color=BLUE, stroke_width=5)
        refl = Arrow(LEFT * 0.8 + DOWN * 0.02, LEFT * 3.2 + UP * 1.25, color=PURPLE, stroke_width=5)
        trans = VGroup(
            Arrow(LEFT * 0.8 + DOWN * 0.2, LEFT * 2.0 + DOWN * 1.35, color=BLUE),
            Arrow(ORIGIN + DOWN * 0.2, RIGHT * 0.3 + DOWN * 1.45, color=TEAL),
            Arrow(RIGHT * 0.7 + DOWN * 0.18, RIGHT * 2.1 + DOWN * 1.25, color=PURPLE),
            Arrow(RIGHT * 1.3 + DOWN * 0.18, RIGHT * 3.2 + DOWN * 1.1, color=ORANGE),
        )
        matrix = formula_card(r"\mathbf{A}[R^E,R^M,T^{Pf},T^{Ps},T^{TM},T^{SV}]^T=\mathbf{b}", "Eq. 46, Appendix B Eq. B1-B7", RED, width=6.3).to_edge(RIGHT, buff=0.6).shift(UP * 0.1)
        bcs = VGroup(txt("位移连续", 20), txt("压力连续", 20), txt("应力条件", 20), txt("切向 H 连续", 20), txt("切向 E 连续", 20)).arrange(DOWN, buff=0.08).next_to(matrix, DOWN, buff=0.18)
        labels = VGroup(txt("fluid side", 20, BLUE).move_to(LEFT * 4.2 + UP * 1.35), txt("porous side", 20, ORANGE).move_to(LEFT * 4.0 + DOWN * 1.45))
        self.play(FadeIn(header), FadeIn(water), FadeIn(porous), ShowCreation(interface), FadeIn(labels), run_time=0.9)
        self.play(GrowArrow(inc), GrowArrow(refl), FadeIn(trans, lag_ratio=0.08), run_time=1.2)
        self.play(FadeIn(matrix), FadeIn(bcs, lag_ratio=0.08), run_time=0.9)
        self.active = VGroup(header, water, porous, interface, inc, refl, trans, labels, matrix, bcs)
        self.hold("interface")

    def normalization_source_segment(self):
        self.clear_active()
        header = self.scene_title("压力归一化和代码声源", "把 2010 系数放进 2011 正演的关键转接")
        cards = VGroup(
            formula_card(r"R^E_{\mathrm{press}}={R^E_{2010}\over\rho_{fl}\omega^2}", "Schakel et al. 2011 Table I / Eq. 3 text", BLUE, width=5.9),
            formula_card(r"\Phi^{TM}_{\mathrm{term}}=-{\alpha^{TM}T^{TM}\over\rho_{fl}\omega^2}", "JAP Table I + code scalar-potential sign", PURPLE, width=5.9),
            formula_card(r"P(\omega)=\int_0^T p_0\,\mathrm{Ricker}(\tau-t_p)e^{-i\omega\tau}d\tau", "code causal pressure source", ORANGE, width=5.9),
            formula_card(r"A(\omega)=P(\omega)R_{\mathrm{ref}}W(f)", "Schakel pressure-source convention + code bandpass", TEAL, width=5.9),
        ).arrange(DOWN, buff=0.18).to_edge(LEFT, buff=0.55).shift(DOWN * 0.05)
        transducer = Circle(radius=0.72).set_fill("#B9D6F2", 1).set_stroke(BLUE, 3).to_edge(RIGHT, buff=2.2).shift(UP * 0.8)
        waves = VGroup(*[Arc(radius=0.9 + 0.33 * i, start_angle=-0.65, angle=1.3, color=BLUE).move_arc_center_to(transducer.get_center()) for i in range(5)])
        dtheta = formula_card(r"D(\theta)={J_1(ka\sin\theta)\over ka\sin\theta}", "Schakel et al. 2011 Eq. 2", RED, width=4.5).next_to(transducer, DOWN, buff=0.65)
        self.play(FadeIn(header), FadeIn(cards, lag_ratio=0.08), run_time=1.4)
        self.play(FadeIn(transducer), ShowCreation(waves, lag_ratio=0.08), FadeIn(dtheta), run_time=1.2)
        self.active = VGroup(header, cards, transducer, waves, dtheta)
        self.hold("normalization_source")

    def sommerfeld_segment(self):
        self.clear_active()
        header = self.scene_title("Schakel 2011 Sommerfeld 积分", "实角分支 + 倏逝分支，把界面系数合成接收点频域电势")
        plane = Axes(x_range=(0, 5, 1), y_range=(0, 4, 1), width=4.5, height=3.2, axis_config={"stroke_color": MUTED, "stroke_width": 2}).to_edge(LEFT, buff=0.75).shift(DOWN * 0.35)
        real = plane.get_graph(lambda x: 0.65 * x + 0.25, x_range=(0.35, 3.5), color=BLUE)
        ev = plane.get_graph(lambda x: 2.6 + 0.10 * np.sin(5 * x), x_range=(0.55, 4.15), color=PURPLE).shift(UP * 0.38)
        labels = VGroup(txt("实角 θ：传播贡献", 21, BLUE).next_to(real, DOWN), txt("γ：倏逝/近场贡献", 21, PURPLE).next_to(ev, UP))
        formula = formula_card(r"\hat{\phi}^{E}=-{iA\over a}\int_0^{\pi/2}J_0J_1e^{ikz_s\cos\theta}R^Ee^{ik_z^Ez_r}d\theta+{A\over a}\int_0^\infty\cdots d\gamma", "Schakel et al. 2011 Geophysics Eq. 5", BLUE, width=7.0, height=1.75, fsize=28).to_edge(RIGHT, buff=0.55).shift(UP * 0.85)
        blocks = VGroup(box("频率 f", TEAL, 1.5), box("角度 θ/γ", BLUE, 1.9), box("R^E 或 TM 项", ORANGE, 2.4), box("接收深度 z", PURPLE, 2.2)).arrange(RIGHT, buff=0.15).next_to(formula, DOWN, buff=0.45)
        wave = self.wave_panel(4.6, 1.55).next_to(blocks, DOWN, buff=0.45)
        self.play(FadeIn(header), FadeIn(plane), ShowCreation(real), ShowCreation(ev), FadeIn(labels), run_time=1.2)
        self.play(FadeIn(formula), FadeIn(blocks, lag_ratio=0.08), run_time=1.1)
        self.play(FadeIn(wave[0]), ShowCreation(wave[1]), ShowCreation(wave[2]), FadeIn(wave[3:]), run_time=1.0)
        self.active = VGroup(header, plane, real, ev, labels, formula, blocks, wave)
        self.hold("sommerfeld")

    def porous_policy_segment(self):
        self.clear_active()
        header = self.scene_title("多孔介质侧：当前代码是单界面模型", "JAP Eq. 8 的有限厚度 slab 被有意识地简化")
        slab = Rectangle(width=4.6, height=3.0).set_fill("#EDE6D6", 1).set_stroke(ORANGE, 3).to_edge(LEFT, buff=0.8).shift(DOWN * 0.15)
        front = Line(slab.get_left() + DOWN * 1.5, slab.get_left() + UP * 1.5, color=RED, stroke_width=5)
        back = DashedLine(slab.get_right() + DOWN * 1.5, slab.get_right() + UP * 1.5, color=MUTED, stroke_width=3)
        tm = Arrow(front.get_center() + LEFT * 0.2, front.get_center() + RIGHT * 2.1, color=PURPLE, stroke_width=5)
        pf = Arrow(front.get_center() + DOWN * 0.45, front.get_center() + RIGHT * 2.0 + DOWN * 0.8, color=BLUE, stroke_width=3)
        cross = VGroup(Line(back.get_center() + UL * 0.35, back.get_center() + DR * 0.35, color=RED), Line(back.get_center() + DL * 0.35, back.get_center() + UR * 0.35, color=RED))
        labels = VGroup(
            mtx(r"\Phi^{TM}_{term}", 30, PURPLE).next_to(tm, UP),
            txt("Pf 共震项：可选，默认关闭", 20, BLUE).next_to(pf, DOWN),
            txt("背界面多次项：当前不实现", 20, MUTED).next_to(back, RIGHT),
        )
        table = VGroup(
            box("JAP Eq. 8：有限厚度实验样品", ORANGE, width=4.4),
            box("当前代码：单界面 RT-SE 正演", TEAL, width=4.4),
            box("目标：突出 interface EM response", RED, width=4.4),
        ).arrange(DOWN, buff=0.22).to_edge(RIGHT, buff=0.75)
        self.play(FadeIn(header), FadeIn(slab), ShowCreation(front), ShowCreation(back), run_time=0.9)
        self.play(GrowArrow(tm), GrowArrow(pf), FadeIn(labels[:2]), run_time=1.0)
        self.play(ShowCreation(cross), FadeIn(labels[2]), FadeIn(table, lag_ratio=0.10), run_time=1.1)
        self.active = VGroup(header, slab, front, back, tm, pf, cross, labels, table)
        self.hold("porous_policy")

    def time_domain_segment(self):
        self.clear_active()
        header = self.scene_title("回到时间域，并明确哪些点参与定量", "正频率反变换、T0、接收侧和 z=0 插值行")
        inv = formula_card(r"u(z,t)=2\mathrm{Re}\int_{f_{min}}^{f_{max}}\hat{u}(z,f)e^{i2\pi ft}df", "positive-frequency inverse transform", BLUE, width=5.6).to_edge(LEFT, buff=0.6).shift(UP * 1.2)
        t0 = formula_card(r"T_0={z_s\over v_f},\quad v_f=\sqrt{K_{fl}/\rho_{fl}}", "interface acoustic arrival", ORANGE, width=5.6).next_to(inv, DOWN, buff=0.3)
        z0 = formula_card(r"u(0,t)={1\over2}\left[u(z_-,t)+u(z_+,t)\right]", "plot-only interpolation; excluded from peaks", RED, width=5.6).next_to(t0, DOWN, buff=0.3)
        wave = self.wave_panel(5.3, 2.25).to_edge(RIGHT, buff=0.65).shift(UP * 0.25)
        checks = VGroup(box("流体侧接收点", BLUE, 2.4), box("多孔侧接收点", RED, 2.4), box("界面插值行排除", ORANGE, 2.4)).arrange(RIGHT, buff=0.22).to_edge(DOWN, buff=0.55)
        self.play(FadeIn(header), FadeIn(inv), FadeIn(t0), FadeIn(z0), run_time=1.2)
        self.play(FadeIn(wave[0]), ShowCreation(wave[1]), ShowCreation(wave[2]), FadeIn(wave[3:]), FadeIn(checks, lag_ratio=0.10), run_time=1.3)
        self.active = VGroup(header, inv, t0, z0, wave, checks)
        self.hold("time_domain")

    def pe_results_segment(self):
        self.clear_active()
        header = self.scene_title("把公式链放回三组 Pe 结果", "连通性增强，但界面 EM 峰值下降")
        axes = Axes(x_range=(0, 16, 4), y_range=(-4.3, 0.3, 1), width=6.2, height=3.45, axis_config={"stroke_color": MUTED, "stroke_width": 2}).to_edge(RIGHT, buff=0.7).shift(DOWN * 0.1)
        ylab = txt("log10 normalized peak", 18, MUTED).next_to(axes, LEFT, buff=0.05).rotate(PI / 2)
        xlab = txt("dissolution time (min)", 18, MUTED).next_to(axes, DOWN, buff=0.08)
        colors = {"Pe=0.1": BLUE, "Pe=1": ORANGE, "Pe=10": RED}
        curves = VGroup()
        legends = VGroup()
        for i, (label, color) in enumerate(colors.items()):
            samples = self.pe_data.get(label, {}).get("samples", [])
            curves.add(self.data_curve(axes, samples, "Amax_waveform_schakel2011_RE_norm", color, ylog=True))
            legends.add(VGroup(Line(ORIGIN, RIGHT * 0.35, color=color, stroke_width=4), txt(label, 18, color, width=1.5)).arrange(RIGHT, buff=0.08))
        legends.arrange(DOWN, buff=0.12).next_to(axes, UP, buff=0.12).align_to(axes, RIGHT)
        stats = VGroup(
            box("孔隙率/渗透率上升", TEAL, width=3.0),
            box("σ(ω) 可上升", PURPLE, width=2.2),
            box("电动耦合下降", BLUE, width=2.4),
            box("界面系数下降", ORANGE, width=2.4),
            box("peak ↓↓↓", RED, width=2.0),
        ).arrange(DOWN, buff=0.16).to_edge(LEFT, buff=0.95).shift(UP * 0.18)
        claim = txt("反直觉结论：水力连通性增强，不必然增强界面电流不平衡。", 21, RED, width=4.25).next_to(stats, DOWN, buff=0.42).align_to(stats, LEFT)
        self.play(FadeIn(header), FadeIn(stats, lag_ratio=0.10), FadeIn(claim), run_time=1.1)
        self.play(FadeIn(axes), FadeIn(ylab), FadeIn(xlab), FadeIn(legends), ShowCreation(curves, lag_ratio=0.10), run_time=1.8)
        self.active = VGroup(header, axes, ylab, xlab, curves, legends, stats, claim)
        self.hold("pe_results")

    def waveform_segment(self):
        self.clear_active()
        header = self.scene_title("波形快照：信号定位在界面附近", "T0 后窄脉冲、上下极性反转、离界面衰减")
        panel = self.wave_panel(5.2, 2.35).to_edge(RIGHT, buff=0.75).shift(UP * 0.75)
        interface = Line(LEFT * 2.4, RIGHT * 2.4, color=INK, stroke_width=3).to_edge(LEFT, buff=0.75).shift(UP * 0.25)
        top = Rectangle(width=4.8, height=1.35).set_fill("#DDEFF8", 1).set_stroke(width=0).next_to(interface, UP, buff=0)
        bottom = Rectangle(width=4.8, height=1.35).set_fill("#EDE6D6", 1).set_stroke(width=0).next_to(interface, DOWN, buff=0)
        dipole = VGroup(Arrow(LEFT * 0.0 + DOWN * 0.6, LEFT * 0.0 + UP * 0.6, color=RED, stroke_width=5), txt("+", 26, RED).move_to(UP * 0.9), txt("-", 26, BLUE).move_to(DOWN * 0.9)).move_to(interface.get_center())
        decay = VGroup(*[Line(interface.get_center() + RIGHT * (0.55 + i * 0.35), interface.get_center() + RIGHT * (0.55 + i * 0.35) + UP * (0.75 - i * 0.12), color=PURPLE, stroke_width=max(2, 6 - i)) for i in range(5)])
        bullets = VGroup(
            box("T0 后出现主脉冲", ORANGE, width=3.0),
            box("z<0 与 z>0 极性相反", PURPLE, width=3.3),
            box("峰值离界面越远越弱", BLUE, width=3.2),
        ).arrange(DOWN, buff=0.18).to_edge(DOWN, buff=0.55)
        self.play(FadeIn(header), FadeIn(top), FadeIn(bottom), ShowCreation(interface), run_time=0.8)
        self.play(FadeIn(dipole), ShowCreation(decay, lag_ratio=0.08), FadeIn(panel[0]), ShowCreation(panel[1]), ShowCreation(panel[2]), FadeIn(panel[3:]), run_time=1.5)
        self.play(FadeIn(bullets, lag_ratio=0.10), run_time=0.8)
        self.active = VGroup(header, top, bottom, interface, dipole, decay, panel, bullets)
        self.hold("waveform")

    def liu_context_segment(self):
        self.clear_active()
        header = self.scene_title("Liu 2018 在这里是什么角色", "有限偏移 VSEP 背景，而不是当前零偏移主积分")
        liu = formula_card(r"u_r,u_t=\iint A(k,\omega)\,R_u/T_u\,e^{-i(\cdots-\omega t)}\,dk\,d\omega", "Liu et al. 2018 Eq. 1-2, finite-offset VSEP", PURPLE, width=6.6).to_edge(LEFT, buff=0.65).shift(UP * 0.55)
        schakel = formula_card(r"\hat{\phi}^{E}=-{iA\over a}\int_\theta(\cdots)d\theta+{A\over a}\int_\gamma(\cdots)d\gamma", "current code: Schakel et al. 2011 Sommerfeld, offset_D=0", BLUE, width=6.6).to_edge(RIGHT, buff=0.65).shift(UP * 0.55)
        bridge = Arrow(liu.get_right(), schakel.get_left(), color=ORANGE)
        note = VGroup(
            box("Liu：解释有限偏移空间分布", PURPLE, width=4.2),
            box("本脚本：解释零偏移界面 EM 随溶蚀衰减", BLUE, width=5.3),
            box("未来扩展：把 Pe 结果放回有限偏移 VSEP", ORANGE, width=5.2),
        ).arrange(DOWN, buff=0.25).to_edge(DOWN, buff=0.8)
        self.play(FadeIn(header), FadeIn(liu), FadeIn(schakel), run_time=1.0)
        self.play(GrowArrow(bridge), FadeIn(note, lag_ratio=0.10), run_time=1.1)
        self.active = VGroup(header, liu, schakel, bridge, note)
        self.hold("liu_context")

    def appendix_material_segment(self):
        self.clear_active()
        header = self.scene_title("公式附录 A：材料参数完整链", "把 L(ω)、电导率拆分和介电混合逐项补齐")
        left = VGroup(
            formula_card(r"L=-{\phi\over\alpha_\infty}{\varepsilon_0\varepsilon_f\zeta\over\eta}\left(1-{2d\over\Lambda}\right)\left[1+{2i\omega\over M\omega_t}\left(1-{2d\over\Lambda}\right)^2\left(1+d\sqrt{{i\omega\rho_f\over\eta}}\right)^2\right]^{-1/2}", "Schakel & Smeulders 2010 Eq. A4", BLUE, width=6.35, height=1.45, fsize=22),
            formula_card(r"\sigma_f=\sum_l(ez_l)^2b_lN_l", "Eq. A7", TEAL, width=6.35, height=0.95, fsize=30),
            formula_card(r"C_{em}=2d\sum_l(ez_l)^2b_lN_l\left[\exp\left(-{ez_l\zeta\over2k_BT}\right)-1\right]", "Eq. A8", ORANGE, width=6.35, height=1.10, fsize=24),
        ).arrange(DOWN, buff=0.13).to_edge(LEFT, buff=0.45).shift(DOWN * 0.10)
        right = VGroup(
            formula_card(r"P_{os}={8k_BTd^2\over\varepsilon_0\varepsilon_f\zeta^2}\sum_lN_l\left[\exp\left(-{ez_l\zeta\over2k_BT}\right)-1\right]", "Eq. A10", PURPLE, width=5.9, height=1.10, fsize=23),
            formula_card(r"C_{os}(\omega)={(\varepsilon_0\varepsilon_f)^2\zeta^2\over2d\eta}P_{os}\left[1+{2d\over P_{os}}\sqrt{{i\omega\rho_f\over\eta}}\right]^{-1}", "Eq. A9", PURPLE, width=5.9, height=1.20, fsize=22),
            formula_card(r"\sigma(\omega)={\phi\sigma_f\over\alpha_\infty}\left[1+{2(C_{em}+C_{os})\over\sigma_f\Lambda}\right]", "Eq. A6", RED, width=5.9, height=0.95, fsize=27),
            formula_card(r"\varepsilon=\varepsilon_0\left[{\phi(\varepsilon_f-\varepsilon_s)\over\alpha_\infty}+\varepsilon_s\right],\quad \bar{\varepsilon}=\varepsilon-i{\sigma\over\omega}+i{\eta L^2\over\omega k}", "Eq. 7 text + Eq. 20", BLUE, width=5.9, height=1.05, fsize=23),
        ).arrange(DOWN, buff=0.12).to_edge(RIGHT, buff=0.45).shift(DOWN * 0.10)
        flow = VGroup(box("ions", TEAL, 1.2), box("ζ,d,Λ", ORANGE, 1.45), box("L,σ", BLUE, 1.45), box("ε̄", RED, 1.1)).arrange(RIGHT, buff=0.14).to_edge(DOWN, buff=0.28)
        arrows = VGroup(*[arrow_between(flow[i], flow[i + 1], MUTED) for i in range(len(flow) - 1)])
        self.play(FadeIn(header), FadeIn(left, lag_ratio=0.06), run_time=1.25)
        self.play(FadeIn(right, lag_ratio=0.06), run_time=1.25)
        self.play(FadeIn(flow, lag_ratio=0.08), ShowCreation(arrows, lag_ratio=0.08), run_time=0.8)
        self.active = VGroup(header, left, right, flow, arrows)
        self.hold("appendix_material")

    def appendix_modes_boundary_segment(self):
        self.clear_active()
        header = self.scene_title("公式附录 B：波模慢度与 B2-B7 矩阵行", "把 rho、E_K、二次根、alpha/beta 和边界矩阵逐项对上")
        density = VGroup(
            formula_card(r"\rho_{12}=\phi\rho_f\left[1+i{\phi\eta\over\omega\rho_fk(\omega)}\right],\quad \rho_{11}=(1-\phi)\rho_s-\rho_{12},\quad \rho_{22}=\phi\rho_f-\rho_{12}", "Eq. 11-13", BLUE, width=6.25, height=1.10, fsize=22),
            formula_card(r"E_K={\eta^2\phi^2L^2\over k^2\bar{\varepsilon}\omega^2},\quad \bar{\rho}_{11}=\rho_{11}-E_K,\ \bar{\rho}_{12}=\rho_{12}+E_K,\ \bar{\rho}_{22}=\rho_{22}-E_K", "Eq. 26-29", TEAL, width=6.25, height=1.10, fsize=22),
            formula_card(r"s_l^2={-d_1\over2d_2}\mp{1\over2}\sqrt{\left({d_1\over d_2}\right)^2-4{d_0\over d_2}}", "P waves Eq. 24-25; TM/SV Eq. 30 with transverse d0-d2", ORANGE, width=6.25, height=1.00, fsize=26),
            formula_card(r"\beta={\mathrm{fluid\ potential}\over\mathrm{solid\ potential}},\quad \alpha={\mathrm{electric\ potential}\over\mathrm{mechanical\ potential}}", "Eq. 36-39", PURPLE, width=6.25, height=0.90, fsize=25),
        ).arrange(DOWN, buff=0.12).to_edge(LEFT, buff=0.42).shift(DOWN * 0.10)
        rows = VGroup(
            formula_card(r"B2:\ [0,k_3^{fl},k_3^{Pf}(1-\phi+\phi\beta^{Pf}),\ldots]", "volume displacement", BLUE, width=5.85, height=0.74, fsize=20, csize=14),
            formula_card(r"B3:\ [0,-\phi\rho_f,(Q+R\beta^{Pf})s_{Pf}^2,\ldots]", "pressure", TEAL, width=5.85, height=0.74, fsize=20, csize=14),
            formula_card(r"B4,B5:\ k_1k_3,\ k_1^2-\omega^2s^2N/(2G),\ -k_1k_3", "stress rows with N1,N2", ORANGE, width=5.85, height=0.82, fsize=19, csize=14),
            formula_card(r"B6:\ [-s_E^2/\mu,0,0,0,\alpha^{TM}s_{TM}^2/\mu,\alpha^{SV}s_{SV}^2/\mu]", "tangential H", PURPLE, width=5.85, height=0.82, fsize=18, csize=14),
            formula_card(r"B7:\ [-k_3^E,0,k_1\alpha^{Pf},k_1\alpha^{Ps},-k_3^{TM}\alpha^{TM},-k_3^{SV}\alpha^{SV}]", "tangential E", RED, width=5.85, height=0.82, fsize=18, csize=14),
            formula_card(r"N_{1,2}=P-{Q(1-\phi)\over\phi}+\left[Q-{R(1-\phi)\over\phi}\right]\beta^{Pf,Ps}", "Appendix B stress combination", INK, width=5.85, height=0.82, fsize=20, csize=14),
        ).arrange(DOWN, buff=0.06).to_edge(RIGHT, buff=0.38).shift(DOWN * 0.08)
        self.play(FadeIn(header), FadeIn(density, lag_ratio=0.06), run_time=1.35)
        self.play(FadeIn(rows, lag_ratio=0.06), run_time=1.5)
        self.active = VGroup(header, density, rows)
        self.hold("appendix_modes_boundary")

    def appendix_source_diagnostics_segment(self):
        self.clear_active()
        header = self.scene_title("公式附录 C：声源、频率网格和诊断", "这些是代码数值实现，不应误写成 Schakel 理论新增公式")
        left = VGroup(
            formula_card(r"p(t)=p_0[1-2(\pi f_0\tau)^2]e^{-(\pi f_0\tau)^2},\quad P(\omega)=\int_0^Tp(\tau)e^{-i\omega\tau}d\tau", "default causal Ricker source", BLUE, width=6.1, height=1.10, fsize=22),
            formula_card(r"A(\omega)=P(\omega)R_{ref}W(f),\quad W(f):\mathrm{half\ cosine\ tapers}", "code bandpass taper", TEAL, width=6.1, height=0.92, fsize=24),
            formula_card(r"R_{ramp}(n)={1\over2}\left[1-\cos{\pi n\over N_{ramp}-1}\right]", "source start ramp", ORANGE, width=6.1, height=0.92, fsize=25),
            formula_card(r"J_0(x)\approx{1\over N}\sum e^{ix\cos\theta_j},\quad J_1(x)\approx{1\over\pi}\int_0^\pi\cos(\theta-x\sin\theta)d\theta", "fallback only if scipy.special unavailable", PURPLE, width=6.1, height=1.02, fsize=21),
        ).arrange(DOWN, buff=0.11).to_edge(LEFT, buff=0.45).shift(DOWN * 0.08)
        right = VGroup(
            formula_card(r"f_j=f_{min}+j{f_{max}-f_{min}\over N_f-1},\quad u(z,t)=2\mathrm{Re}\int\hat{u}(z,f)e^{i2\pi ft}df", "linear grid + positive-frequency inverse transform", BLUE, width=5.95, height=1.03, fsize=22),
            formula_card(r"X_{norm}(t_d)={X(t_d)\over X(t_{d,0})},\quad A_{max}=\max_{z\ne0,t}|u(z,t)|", "first valid normalization + peak", RED, width=5.95, height=0.96, fsize=23),
            formula_card(r"\Delta A_{rel}^{(N)}={A^{(N)}-A^{(N_{prev})}\over A^{(N_{prev})}},\quad \Delta t_{peak}=t_{peak}^{(N)}-t_{peak}^{(N_{prev})}", "convergence diagnostics", ORANGE, width=5.95, height=1.00, fsize=21),
            formula_card(r"u(-d,t^*)u(+d,t^*)<0,\quad { \max_{t<T_0}|u| \over \max_{t\ge T_0}|u| }", "polarity reversal + pre/post T0 ratio", PURPLE, width=5.95, height=0.96, fsize=23),
        ).arrange(DOWN, buff=0.12).to_edge(RIGHT, buff=0.45).shift(DOWN * 0.08)
        note = box("fig4_digitized：Schakel 2011 Fig. 4 的视觉近似；默认仍是 Ricker", RED, width=8.8, height=0.55, size=18).to_edge(DOWN, buff=0.23)
        self.play(FadeIn(header), FadeIn(left, lag_ratio=0.06), run_time=1.35)
        self.play(FadeIn(right, lag_ratio=0.06), FadeIn(note), run_time=1.45)
        self.active = VGroup(header, left, right, note)
        self.hold("appendix_source_diagnostics")

    def ending_segment(self):
        self.clear_active()
        header = self.scene_title("最终机制链", "从反应输运到 JGR Solid Earth 可讲的震电故事")
        steps = VGroup(
            box("溶蚀", TEAL, 1.35),
            math_box(r"\phi,k_0,\alpha_\infty,c_H", TEAL, 2.2, size=24),
            math_box(r"k,L,\sigma,\bar{\varepsilon}", BLUE, 2.0, size=24),
            math_box(r"R^E,T^{TM}", ORANGE, 2.0, size=24),
            box("Sommerfeld", PURPLE, 2.2),
            box("peak ↓", RED, 1.65),
        ).arrange(RIGHT, buff=0.16).move_to(UP * 1.25)
        arrows = VGroup(*[arrow_between(steps[i], steps[i + 1], color=MUTED) for i in range(len(steps) - 1)])
        formula = mtx(r"\{\phi,k_0,\alpha_\infty,c_H\}(t_d)\rightarrow\{L,\sigma,\bar{\varepsilon}\}\rightarrow R^E,T^{TM}\rightarrow u(z,t)\rightarrow A_{\max}(t_d)", 34, INK, width=11.2).next_to(steps, DOWN, buff=0.75)
        takeaways = VGroup(
            txt("Schakel 2010：材料和界面转换系数。", 24, width=9.5),
            txt("Schakel 2011：压力源和 Sommerfeld 波形合成。", 24, width=9.5),
            txt("你的结果：Pe 越高，峰值越快跨数量级降低。", 24, RED, width=9.5),
        ).arrange(DOWN, buff=0.18).next_to(formula, DOWN, buff=0.55)
        self.play(FadeIn(header), FadeIn(steps, lag_ratio=0.08), run_time=1.0)
        self.play(ShowCreation(arrows, lag_ratio=0.10), FadeIn(formula), run_time=1.2)
        self.play(FadeIn(takeaways, lag_ratio=0.12), run_time=1.0)
        self.active = VGroup(header, steps, arrows, formula, takeaways)
        self.hold("ending")
