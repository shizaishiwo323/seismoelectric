from __future__ import annotations

import json
from pathlib import Path

from manimlib import *


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "hyperframes" / "schakel2011_sommerfeld_explainer" / "assets"
TIMELINE_PATH = ASSETS / "timeline.json"

FONT = "Microsoft YaHei"
INK = "#263238"
MUTED = "#607D8B"
BLUE = "#2E86AB"
TEAL = "#0EAD69"
ORANGE = "#F28C28"
RED = "#C44536"
PURPLE = "#7B61FF"
YELLOW = "#F7C948"
PAPER = "#F7F9FC"


def fit(mob, width=12.0, height=6.3):
    if mob.get_width() > width:
        mob.set_width(width)
    if mob.get_height() > height:
        mob.set_height(height)
    return mob


def txt(s, size=30, color=INK, width=12.0):
    mob = Text(s, font=FONT, font_size=size)
    mob.set_color(color)
    return fit(mob, width=width)


def mtx(s, size=42, color=INK, width=12.0):
    mob = Tex(s, font_size=size)
    mob.set_color(color)
    return fit(mob, width=width)


def box(label, color=BLUE, width=2.3, height=0.72, size=24):
    rect = RoundedRectangle(width=width, height=height, corner_radius=0.08)
    rect.set_fill(color, opacity=0.16).set_stroke(color, width=2)
    label_mob = txt(label, size=size, color=INK, width=width - 0.18)
    return VGroup(rect, label_mob)


def formula_card(formula, caption, color=BLUE):
    card = RoundedRectangle(width=5.3, height=1.35, corner_radius=0.1)
    card.set_fill(WHITE, opacity=0.96).set_stroke(color, width=2)
    f = mtx(formula, size=34, color=INK, width=4.95)
    c = txt(caption, size=20, color=MUTED, width=4.95)
    group = VGroup(f, c).arrange(DOWN, buff=0.12)
    return VGroup(card, group)


def arrow_between(left, right, color=MUTED):
    return Arrow(left.get_right() + 0.08 * RIGHT, right.get_left() + 0.08 * LEFT, buff=0.06, color=color)


class Schakel2011SommerfeldExplainer(Scene):
    CONFIG = {
        "camera_config": {"background_color": PAPER, "pixel_width": 1280, "pixel_height": 720},
    }

    def setup(self):
        if TIMELINE_PATH.exists():
            data = json.loads(TIMELINE_PATH.read_text(encoding="utf-8"))
            self.timeline = {item["key"]: item for item in data["segments"]}
            self.total_duration = data["total_duration"]
        else:
            keys = [
                "title",
                "two_times",
                "rt_inputs",
                "dynamic_bridge",
                "elastic_modes",
                "interface_matrix",
                "coefficients",
                "source",
                "sommerfeld",
                "porous_side",
                "time_domain",
                "summary",
            ]
            self.timeline = {key: {"end": (i + 1) * 22.0} for i, key in enumerate(keys)}
            self.total_duration = len(keys) * 22.0
        self.active = VGroup()

    def hold(self, key):
        target = float(self.timeline[key]["end"])
        remaining = max(0.2, target - self.time)
        self.wait(remaining)

    def clear_active(self):
        if len(self.active) > 0:
            self.play(FadeOut(self.active), run_time=0.55)
        self.active = VGroup()

    def scene_title(self, title, subtitle=None):
        t = txt(title, size=34, color=INK, width=11.5).to_edge(UP, buff=0.32)
        underline = Line(LEFT * 5.6, RIGHT * 5.6, color=BLUE).next_to(t, DOWN, buff=0.12)
        group = VGroup(t, underline)
        if subtitle:
            sub = txt(subtitle, size=20, color=MUTED, width=11.5).next_to(underline, DOWN, buff=0.12)
            group.add(sub)
        self.play(FadeIn(t, shift=0.15 * DOWN), ShowCreation(underline), run_time=0.8)
        if subtitle:
            self.play(FadeIn(sub), run_time=0.4)
        return group

    def vector_micro_pores(self):
        matrix = RoundedRectangle(width=4.8, height=3.15, corner_radius=0.08)
        matrix.set_fill("#E6EEF5", opacity=1).set_stroke("#9AB4C8", width=2)
        pores = VGroup()
        centers = [(-1.8, 0.75), (-0.7, 0.2), (0.65, 0.75), (1.65, -0.1), (-1.4, -0.85), (0.2, -0.78)]
        radii = [0.38, 0.48, 0.35, 0.45, 0.36, 0.42]
        for (x, y), r in zip(centers, radii):
            c = Circle(radius=r).set_fill(WHITE, opacity=1).set_stroke(BLUE, width=2)
            c.move_to(matrix.get_center() + x * RIGHT + y * UP)
            pores.add(c)
        channels = VGroup()
        for a, b in [(0, 1), (1, 2), (2, 3), (1, 5), (4, 5)]:
            line = Line(pores[a].get_center(), pores[b].get_center(), color=BLUE, stroke_width=10)
            line.set_opacity(0.32)
            channels.add(line)
        ions = VGroup()
        for x, y, sign, color in [
            (-1.8, 0.75, "+", RED),
            (-0.7, 0.2, "-", TEAL),
            (0.65, 0.75, "+", RED),
            (1.65, -0.1, "-", TEAL),
            (0.2, -0.78, "+", RED),
        ]:
            ions.add(txt(sign, size=26, color=color).move_to(matrix.get_center() + x * RIGHT + y * UP))
        label = txt("溶蚀改变孔隙连通性和流体化学", size=22, color=INK, width=4.6).next_to(matrix, DOWN, buff=0.18)
        return VGroup(matrix, channels, pores, ions, label)

    def wave_panel(self, width=5.0, height=2.2):
        axes = Axes(
            x_range=(0, 6, 1),
            y_range=(-1.2, 1.2, 1),
            width=width,
            height=height,
            axis_config={"stroke_color": MUTED, "stroke_width": 2},
        )
        curve1 = axes.get_graph(lambda x: 0.78 * np.exp(-0.38 * (x - 3.1) ** 2) * np.sin(5.3 * x), color=BLUE)
        curve2 = axes.get_graph(lambda x: -0.58 * np.exp(-0.55 * (x - 3.2) ** 2) * np.sin(5.3 * x), color=RED)
        t0 = DashedLine(axes.c2p(2.55, -1.05), axes.c2p(2.55, 1.05), color=ORANGE)
        t0_label = VGroup(mtx(r"T_0", size=28, color=ORANGE), txt("界面到时", size=18, color=ORANGE)).arrange(RIGHT, buff=0.08)
        t0_label.next_to(t0, UP, buff=0.05)
        return VGroup(axes, curve1, curve2, t0, t0_label)

    def construct(self):
        bg = Rectangle(width=FRAME_WIDTH + 1, height=FRAME_HEIGHT + 1)
        bg.set_fill(PAPER, opacity=1).set_stroke(width=0)
        self.add(bg)
        self.title_segment()
        self.two_times_segment()
        self.rt_inputs_segment()
        self.dynamic_bridge_segment()
        self.elastic_modes_segment()
        self.interface_matrix_segment()
        self.coefficients_segment()
        self.source_segment()
        self.sommerfeld_segment()
        self.porous_side_segment()
        self.time_domain_segment()
        self.summary_segment()
        self.wait(max(0.2, self.total_duration - self.time))

    def title_segment(self):
        self.clear_active()
        title = txt("Schakel 2011 Sommerfeld 正演公式讲解", size=42, color=INK, width=11)
        formula = mtx(r"\{\phi,k_0,\alpha_\infty,c_H\}\rightarrow \hat{\phi}^{E}(\omega,z)\rightarrow u(z,t)", size=42, color=BLUE)
        subtitle = txt("从反应输运参数到 interface EM response 波形", size=26, color=MUTED, width=10.5)
        chain = VGroup(title, formula, subtitle).arrange(DOWN, buff=0.32).move_to(ORIGIN + 0.35 * UP)
        left = self.vector_micro_pores().scale(0.62).to_edge(LEFT, buff=0.55).shift(0.8 * DOWN)
        wave = self.wave_panel(width=4.3, height=1.85).to_edge(RIGHT, buff=0.65).shift(1.15 * DOWN)
        arrow = Arrow(left.get_right(), wave.get_left(), color=ORANGE, stroke_width=6)
        self.play(FadeIn(chain, shift=0.25 * UP), run_time=1.0)
        self.play(FadeIn(left), GrowArrow(arrow), ShowCreation(wave[1]), ShowCreation(wave[2]), FadeIn(wave[0]), run_time=1.6)
        self.active = VGroup(chain, left, arrow, wave)
        self.hold("title")

    def two_times_segment(self):
        self.clear_active()
        header = self.scene_title("先分清两条时间轴", "慢变量更新材料，快变量合成波形")
        slow = NumberLine(x_range=(0, 10, 2), width=5.3, color=TEAL).shift(UP * 1.1 + LEFT * 2.2)
        fast = NumberLine(x_range=(0, 10, 2), width=5.3, color=BLUE).shift(DOWN * 1.2 + RIGHT * 2.2)
        slow_label = txt("dissolution time / Time_s：秒到小时", size=28, color=TEAL, width=5.5).next_to(slow, UP)
        fast_label = txt("waveform time：微秒级传播时间", size=28, color=BLUE, width=5.4).next_to(fast, DOWN)
        state_box = box("选定一个溶蚀状态", TEAL, width=3.1).next_to(slow, DOWN, buff=0.3)
        response_box = box("在该状态下算频域响应", BLUE, width=3.35).next_to(fast, UP, buff=0.3)
        connector = Arrow(state_box.get_right(), response_box.get_left(), color=ORANGE, buff=0.08)
        td = mtx(r"t_d", size=34, color=TEAL).next_to(slow, RIGHT)
        tw = mtx(r"t", size=34, color=BLUE).next_to(fast, RIGHT)
        phase = VGroup(
            formula_card(r"e^{i\omega t},\quad e^{-i\mathbf{k}\cdot\mathbf{x}}", "相位约定", ORANGE),
            formula_card(r"z<0:\mathrm{fluid},\quad z>0:\mathrm{porous}", "坐标和衰减分支", PURPLE),
        ).arrange(RIGHT, buff=0.22).scale(0.68).to_edge(DOWN, buff=0.28)
        self.play(FadeIn(header), ShowCreation(slow), ShowCreation(fast), FadeIn(slow_label), FadeIn(fast_label))
        self.play(FadeIn(state_box), GrowArrow(connector), FadeIn(response_box), FadeIn(td), FadeIn(tw), run_time=1.4)
        self.play(FadeIn(phase, lag_ratio=0.12), run_time=0.8)
        self.active = VGroup(header, slow, fast, slow_label, fast_label, state_box, response_box, connector, td, tw, phase)
        self.hold("two_times")

    def rt_inputs_segment(self):
        self.clear_active()
        header = self.scene_title("反应输运输出怎样进入模型", "单位换算和 zeta 映射是 RT-SE 接口假设")
        inputs = VGroup(
            box("孔隙率", TEAL),
            box("渗透率", TEAL),
            box("曲折度", TEAL),
            box("H+ 浓度", TEAL),
        ).arrange(DOWN, buff=0.18).to_edge(LEFT, buff=0.72).shift(DOWN * 0.2)
        formulas = VGroup(
            formula_card(r"k_0[\mathrm{m^2}]=k_0[\mathrm{mD}]\times 9.869233\times10^{-16}", "渗透率单位换算", BLUE),
            formula_card(r"pH=-\log_{10}\max(c_H,10^{-7})", "H+ 浓度下限与 pH", BLUE),
            formula_card(r"\zeta=\left[0.010+0.025\log_{10}C\right]\frac{pH-2}{5}", "Schakel 2010 Eq. A5", ORANGE),
        ).arrange(DOWN, buff=0.22).to_edge(RIGHT, buff=0.58).shift(DOWN * 0.08)
        arrows = VGroup(*[Arrow(inputs.get_right(), formulas[i].get_left(), color=MUTED, buff=0.15) for i in range(3)])
        note = txt("C 的取法和上覆流体电导率开关，是项目级建模选择", size=24, color=RED, width=6.3).next_to(formulas, DOWN, buff=0.22)
        self.play(FadeIn(header), FadeIn(inputs, lag_ratio=0.12), run_time=1.0)
        self.play(GrowArrow(arrows[0]), FadeIn(formulas[0]), GrowArrow(arrows[1]), FadeIn(formulas[1]), run_time=1.3)
        self.play(GrowArrow(arrows[2]), FadeIn(formulas[2]), FadeIn(note), run_time=1.1)
        self.active = VGroup(header, inputs, formulas, arrows, note)
        self.hold("rt_inputs")

    def dynamic_bridge_segment(self):
        self.clear_active()
        header = self.scene_title("动态渗透率、耦合系数和电导率", "这是溶蚀影响界面 EM response 的第一座桥")
        pore = self.vector_micro_pores().scale(0.8).to_edge(LEFT, buff=0.65).shift(DOWN * 0.35)
        cards = VGroup(
            formula_card(r"k(\omega)=k_0\left[\sqrt{1+i{\omega M\over2\omega_t}}+i{\omega\over\omega_t}\right]^{-1}", "动态渗透率", BLUE),
            formula_card(r"L(\omega)\propto -{\phi\over\alpha_\infty}{\varepsilon_0\varepsilon_f\zeta\over\eta}", "电动耦合系数", TEAL),
            formula_card(r"\bar{\varepsilon}=\varepsilon-i{\sigma(\omega)\over\omega}+i{\eta L^2\over\omega k}", "有效介电常数", PURPLE),
        ).arrange(DOWN, buff=0.22).to_edge(RIGHT, buff=0.58).shift(DOWN * 0.15)
        flow = VGroup(
            Arrow(pore.get_right(), cards[0].get_left(), color=BLUE),
            Arrow(cards[0].get_bottom(), cards[1].get_top(), color=TEAL),
            Arrow(cards[1].get_bottom(), cards[2].get_top(), color=PURPLE),
        )
        sigma_note = formula_card(r"\sigma(\omega)\ne\sigma_{fl}", "多孔体动态电导率 vs 上覆自由流体电导率", RED)
        sigma_note.scale(0.72).to_edge(DOWN, buff=0.28).shift(RIGHT * 2.35)
        edl = Circle(radius=0.85, color=ORANGE).move_to(pore[2][1].get_center())
        edl.set_stroke(ORANGE, width=4, opacity=0.75)
        self.play(FadeIn(header), FadeIn(pore), run_time=1.0)
        self.play(ShowCreation(edl), GrowArrow(flow[0]), FadeIn(cards[0]), run_time=1.1)
        self.play(GrowArrow(flow[1]), FadeIn(cards[1]), Rotate(edl, angle=TAU, run_time=1.0), run_time=1.2)
        self.play(GrowArrow(flow[2]), FadeIn(cards[2]), run_time=1.0)
        self.play(FadeIn(sigma_note), run_time=0.6)
        self.active = VGroup(header, pore, cards, flow, edl, sigma_note)
        self.hold("dynamic_bridge")

    def elastic_modes_segment(self):
        self.clear_active()
        header = self.scene_title("孔弹波模：机械波和电磁势怎样绑定", "Biot 系数、有效密度、慢度根、alpha/beta 幅值比")
        left = VGroup(
            formula_card(r"A,Q,R,\quad P=A+2G", "Biot 弹性系数", TEAL),
            formula_card(r"\bar{\rho}_{ij}=\rho_{ij}\pm E_K", "电动耦合修正密度", ORANGE),
        ).arrange(DOWN, buff=0.35).to_edge(LEFT, buff=0.65)
        roots = VGroup(
            box("快 P", BLUE, width=1.4),
            box("慢 P", BLUE, width=1.4),
            box("TM", PURPLE, width=1.4),
            box("SV", PURPLE, width=1.4),
        )
        roots[0].move_to(LEFT * 0.85 + UP * 0.25)
        roots[1].move_to(RIGHT * 0.85 + UP * 0.25)
        roots[2].move_to(LEFT * 0.85 + DOWN * 0.65)
        roots[3].move_to(RIGHT * 0.85 + DOWN * 0.65)
        roots.move_to(ORIGIN + 0.15 * DOWN)
        right = VGroup(
            formula_card(r"\beta={\mathrm{fluid\ potential}\over\mathrm{solid\ potential}}", "流固幅值比", BLUE),
            formula_card(r"\alpha={\mathrm{electric\ potential}\over\mathrm{mechanical\ potential}}", "电场幅值比", PURPLE),
        ).arrange(DOWN, buff=0.35).to_edge(RIGHT, buff=0.65)
        arrows = VGroup(
            Arrow(left.get_right(), roots.get_left(), color=MUTED),
            Arrow(roots.get_right(), right.get_left(), color=MUTED),
        )
        self.play(FadeIn(header), FadeIn(left, lag_ratio=0.15), run_time=1.0)
        self.play(GrowArrow(arrows[0]), FadeIn(roots, lag_ratio=0.15), run_time=1.2)
        self.play(GrowArrow(arrows[1]), FadeIn(right, lag_ratio=0.15), run_time=1.0)
        self.active = VGroup(header, left, roots, right, arrows)
        self.hold("elastic_modes")

    def interface_matrix_segment(self):
        self.clear_active()
        header = self.scene_title("界面边界条件：六个未知量，六个条件", "open-pore 界面把机械场和电磁场锁在一起")
        water = Rectangle(width=10.9, height=2.05).set_fill("#DDEFF8", 1).set_stroke(width=0).to_edge(UP, buff=1.15)
        porous = Rectangle(width=10.9, height=2.05).set_fill("#EDE6D6", 1).set_stroke(width=0).to_edge(DOWN, buff=0.35)
        interface = Line(LEFT * 5.45, RIGHT * 5.45, color=INK, stroke_width=4).move_to(ORIGIN + DOWN * 0.15)
        inc = Arrow(LEFT * 3.9 + UP * 1.4, LEFT * 1.4 + DOWN * 0.08, color=BLUE, stroke_width=5)
        refl_e = Arrow(LEFT * 0.9 + DOWN * 0.02, LEFT * 3.2 + UP * 1.35, color=PURPLE, stroke_width=5)
        trans = VGroup(
            Arrow(LEFT * 0.8 + DOWN * 0.2, LEFT * 2.0 + DOWN * 1.4, color=BLUE),
            Arrow(ORIGIN + DOWN * 0.2, RIGHT * 0.3 + DOWN * 1.55, color=TEAL),
            Arrow(RIGHT * 0.75 + DOWN * 0.18, RIGHT * 2.2 + DOWN * 1.32, color=PURPLE),
            Arrow(RIGHT * 1.35 + DOWN * 0.18, RIGHT * 3.4 + DOWN * 1.2, color=ORANGE),
        )
        labels = VGroup(
            txt("流体侧", size=25, color=BLUE).move_to(LEFT * 4.5 + UP * 1.45),
            txt("多孔介质侧", size=25, color=ORANGE).move_to(LEFT * 4.25 + DOWN * 1.55),
            mtx(r"R^E,\ R^M,\ T^{Pf},T^{Ps},T^{TM},T^{SV}", size=32, color=INK).next_to(interface, UP, buff=0.35),
        )
        matrix = mtx(r"\mathbf{A}\,[R^E,R^M,T^{Pf},T^{Ps},T^{TM},T^{SV}]^T=\mathbf{b}", size=40, color=INK)
        matrix.to_edge(RIGHT, buff=0.62).shift(DOWN * 0.2)
        bc = VGroup(
            txt("位移连续", size=22, color=INK),
            txt("压力连续", size=22, color=INK),
            txt("应力条件", size=22, color=INK),
            txt("切向 H 连续", size=22, color=INK),
            txt("切向 E 连续", size=22, color=INK),
        ).arrange(DOWN, buff=0.12).next_to(matrix, DOWN, buff=0.18)
        self.play(FadeIn(header), FadeIn(water), FadeIn(porous), ShowCreation(interface))
        self.play(GrowArrow(inc), GrowArrow(refl_e), FadeIn(trans, lag_ratio=0.1), FadeIn(labels), run_time=1.6)
        self.play(FadeIn(matrix), FadeIn(bc, lag_ratio=0.08), run_time=1.0)
        self.active = VGroup(header, water, porous, interface, inc, refl_e, trans, labels, matrix, bc)
        self.hold("interface_matrix")

    def coefficients_segment(self):
        self.clear_active()
        header = self.scene_title("位移归一化系数怎样变成压力归一化系数", "2010 的界面矩阵系数进入 2011 Sommerfeld 源项")
        left = formula_card(r"R^E_{2010},\quad T^{TM}_{2010},\quad T^{Pf}_{2010}", "六阶系统直接解出的系数", BLUE).to_edge(LEFT, buff=0.85)
        mid = box("压力归一化", ORANGE, width=2.6, height=0.9, size=26)
        right = VGroup(
            formula_card(r"R^E_{\mathrm{press}}={R^E_{2010}\over\rho_{fl}\omega^2}", "流体侧反射 EM", PURPLE),
            formula_card(r"\Phi^{TM}_{\mathrm{term}}=-{\alpha^{TM}T^{TM}\over\rho_{fl}\omega^2}", "多孔侧 TM 电势项", TEAL),
        ).arrange(DOWN, buff=0.35).to_edge(RIGHT, buff=0.7)
        arrows = VGroup(Arrow(left.get_right(), mid.get_left(), color=MUTED), Arrow(mid.get_right(), right.get_left(), color=MUTED))
        note = txt("前导负号是当前代码的标量电势符号约定；归一化来源是 Schakel 2011 JAP Table I", size=22, color=RED, width=10.5).to_edge(DOWN, buff=0.4)
        self.play(FadeIn(header), FadeIn(left), run_time=0.9)
        self.play(GrowArrow(arrows[0]), FadeIn(mid), GrowArrow(arrows[1]), FadeIn(right, lag_ratio=0.14), run_time=1.4)
        self.play(FadeIn(note), run_time=0.5)
        self.active = VGroup(header, left, mid, right, arrows, note)
        self.hold("coefficients")

    def source_segment(self):
        self.clear_active()
        header = self.scene_title("声源谱与圆形活塞方向性", "实验压力记录可以换成代码中的因果 Ricker 压力源")
        transducer = Circle(radius=0.75).set_fill("#B9D6F2", 1).set_stroke(BLUE, 3).to_edge(LEFT, buff=1.1).shift(UP * 0.15)
        waves = VGroup(*[Arc(radius=1.0 + 0.32 * i, start_angle=-0.7, angle=1.4, color=BLUE).move_arc_center_to(transducer.get_center()) for i in range(5)])
        source_formula = formula_card(r"p(t)=p_0\,\mathrm{Ricker}(t-t_p;f_0)", "代码默认压力源", BLUE).move_to(UP * 1.35 + RIGHT * 1.65)
        spectrum = formula_card(r"A(\omega)=P(\omega)R_{\mathrm{ref}}W(f)", "因果傅里叶谱与带通窗", ORANGE).next_to(source_formula, DOWN, buff=0.32)
        directivity = formula_card(r"D(\theta)={J_1(ka\sin\theta)\over ka\sin\theta}", "圆形活塞方向性", PURPLE).next_to(spectrum, DOWN, buff=0.32)
        theta = Arc(radius=1.05, start_angle=0, angle=0.7, color=ORANGE).move_arc_center_to(transducer.get_center())
        theta_label = mtx(r"\theta", size=30, color=ORANGE).next_to(theta, RIGHT, buff=0.05)
        self.play(FadeIn(header), FadeIn(transducer), ShowCreation(waves, lag_ratio=0.12), run_time=1.2)
        self.play(FadeIn(source_formula), FadeIn(spectrum), FadeIn(directivity), ShowCreation(theta), FadeIn(theta_label), run_time=1.3)
        self.active = VGroup(header, transducer, waves, source_formula, spectrum, directivity, theta, theta_label)
        self.hold("source")

    def sommerfeld_segment(self):
        self.clear_active()
        header = self.scene_title("Sommerfeld 积分：把所有角度和频率叠加起来", "实角分支给传播波，倏逝分支给近场衰减贡献")
        plane = Axes(
            x_range=(0, 5, 1),
            y_range=(0, 4, 1),
            width=4.8,
            height=3.5,
            axis_config={"stroke_color": MUTED, "stroke_width": 2},
        ).to_edge(LEFT, buff=0.75).shift(DOWN * 0.25)
        real_path = plane.get_graph(lambda x: 0.65 * x + 0.25, x_range=(0.35, 3.6), color=BLUE)
        ev_path = plane.get_graph(lambda x: 2.5 + 0.12 * np.sin(5 * x), x_range=(0.6, 4.2), color=PURPLE)
        ev_path.shift(UP * 0.4)
        labels = VGroup(
            txt("实角分支", size=22, color=BLUE).next_to(real_path, DOWN),
            txt("倏逝分支", size=22, color=PURPLE).next_to(ev_path, UP),
            mtx(r"\theta", size=28, color=BLUE).next_to(plane, DOWN),
            mtx(r"\gamma", size=28, color=PURPLE).next_to(plane, RIGHT),
        )
        formula = mtx(
            r"\hat{\phi}^{E}=-{iA\over a}\int_{\theta}J_0J_1e^{ikz_s\cos\theta}R^E e^{ik_z^Ez_r}d\theta+{A\over a}\int_{\gamma}\cdots d\gamma",
            size=32,
            color=INK,
            width=6.7,
        ).to_edge(RIGHT, buff=0.55).shift(UP * 0.85)
        stack = VGroup(
            box("频率", TEAL, width=1.5),
            box("角度", BLUE, width=1.5),
            box("转换系数", ORANGE, width=2.0),
            box("接收深度", PURPLE, width=2.0),
        ).arrange(RIGHT, buff=0.18).next_to(formula, DOWN, buff=0.45)
        result = self.wave_panel(width=4.7, height=1.65).next_to(stack, DOWN, buff=0.42)
        self.play(FadeIn(header), FadeIn(plane), ShowCreation(real_path), ShowCreation(ev_path), FadeIn(labels), run_time=1.3)
        self.play(FadeIn(formula), FadeIn(stack, lag_ratio=0.1), run_time=1.0)
        self.play(ShowCreation(result[1]), ShowCreation(result[2]), FadeIn(result[0]), FadeIn(result[3:]), run_time=1.2)
        self.active = VGroup(header, plane, real_path, ev_path, labels, formula, stack, result)
        self.hold("sommerfeld")

    def porous_side_segment(self):
        self.clear_active()
        header = self.scene_title("多孔介质侧：当前代码是单界面模型", "保留前界面 TM 项，默认不做有限厚度 slab 多次反射")
        slab = Rectangle(width=4.6, height=3.0).set_fill("#EDE6D6", 1).set_stroke(ORANGE, 3).to_edge(LEFT, buff=0.85).shift(DOWN * 0.1)
        front = Line(slab.get_left() + DOWN * 1.5, slab.get_left() + UP * 1.5, color=RED, stroke_width=5)
        back = DashedLine(slab.get_right() + DOWN * 1.5, slab.get_right() + UP * 1.5, color=MUTED, stroke_width=3)
        tm_arrow = Arrow(front.get_center() + LEFT * 0.25, front.get_center() + RIGHT * 2.0, color=PURPLE, stroke_width=5)
        pf_arrow = Arrow(front.get_center() + DOWN * 0.45, front.get_center() + RIGHT * 2.0 + DOWN * 0.85, color=BLUE, stroke_width=3)
        crossed = VGroup(Line(back.get_center() + UL * 0.35, back.get_center() + DR * 0.35, color=RED), Line(back.get_center() + DL * 0.35, back.get_center() + UR * 0.35, color=RED))
        labels = VGroup(
            txt("前界面", size=22, color=RED).next_to(front, LEFT),
            txt("背界面多次项：当前不实现", size=22, color=MUTED).next_to(back, RIGHT),
            mtx(r"\Phi^{TM}_{\mathrm{term}}", size=31, color=PURPLE).next_to(tm_arrow, UP),
            txt("Pf 共震项：诊断开关，默认关闭", size=21, color=BLUE).next_to(pf_arrow, DOWN),
        )
        table = VGroup(
            box("JAP Eq. 8：有限厚度实验样品", ORANGE, width=4.2),
            box("当前代码：单界面 RT-SE 正演", TEAL, width=4.2),
            box("目标：突出 interface EM response", BLUE, width=4.2),
        ).arrange(DOWN, buff=0.22).to_edge(RIGHT, buff=0.7)
        self.play(FadeIn(header), FadeIn(slab), ShowCreation(front), ShowCreation(back), run_time=0.9)
        self.play(GrowArrow(tm_arrow), FadeIn(labels[2]), GrowArrow(pf_arrow), FadeIn(labels[3]), run_time=1.2)
        self.play(FadeIn(labels[0]), FadeIn(labels[1]), ShowCreation(crossed), FadeIn(table, lag_ratio=0.12), run_time=1.2)
        self.active = VGroup(header, slab, front, back, tm_arrow, pf_arrow, crossed, labels, table)
        self.hold("porous_side")

    def time_domain_segment(self):
        self.clear_active()
        header = self.scene_title("从频率响应回到时间域，并检查结果可信度", "正频率积分、T0 标记、峰值和极性诊断")
        inverse = formula_card(r"u(z,t)=2\,\mathrm{Re}\int_{f_{\min}}^{f_{\max}}\hat{u}(z,f)e^{i2\pi ft}\,df", "Schakel 相位约定下的反变换", BLUE)
        inverse.to_edge(LEFT, buff=0.65).shift(UP * 1.05)
        t0 = formula_card(r"T_0={z_s\over v_f}", "界面声波到时", ORANGE).next_to(inverse, DOWN, buff=0.35)
        wave = self.wave_panel(width=5.6, height=2.35).to_edge(RIGHT, buff=0.55).shift(UP * 0.2)
        checks = VGroup(
            box("收敛诊断", TEAL, width=2.1),
            box("上下极性反转", PURPLE, width=2.5),
            box("T0 前后能量比", RED, width=2.6),
        ).arrange(RIGHT, buff=0.25).to_edge(DOWN, buff=0.45)
        interface_note = formula_card(r"z=0\ \mathrm{row}:\ \mathrm{plot\ only}", "界面行仅用于绘图插值，峰值诊断排除", RED)
        interface_note.scale(0.72).next_to(checks, UP, buff=0.25)
        arrows = VGroup(
            Arrow(inverse.get_right(), wave.get_left(), color=MUTED),
            Arrow(t0.get_right(), wave.get_left() + DOWN * 0.45, color=ORANGE),
        )
        self.play(FadeIn(header), FadeIn(inverse), FadeIn(t0), run_time=0.9)
        self.play(GrowArrow(arrows[0]), GrowArrow(arrows[1]), FadeIn(wave[0]), ShowCreation(wave[1]), ShowCreation(wave[2]), FadeIn(wave[3:]), run_time=1.4)
        self.play(FadeIn(interface_note), FadeIn(checks, lag_ratio=0.12), run_time=0.8)
        self.active = VGroup(header, inverse, t0, wave, checks, arrows, interface_note)
        self.hold("time_domain")

    def summary_segment(self):
        self.clear_active()
        header = self.scene_title("最后的论文机制链", "把公式列表变成一条可讲的证据链")
        steps = VGroup(
            box("溶蚀", TEAL, width=1.55),
            box("动态系数", BLUE, width=2.15),
            box("界面矩阵", ORANGE, width=2.05),
            box("压力归一化", PURPLE, width=2.35),
            box("Sommerfeld 积分", RED, width=2.65),
            box("interface EM 波形", TEAL, width=2.75),
        ).arrange(RIGHT, buff=0.18).move_to(UP * 1.1)
        arrows = VGroup(*[arrow_between(steps[i], steps[i + 1], color=MUTED) for i in range(len(steps) - 1)])
        key_formula = mtx(
            r"\{\phi,k_0,\alpha_\infty,c_H\}(t_d)\rightarrow\{k,L,\sigma,\bar{\varepsilon}\}\rightarrow R^E,T^{TM}\rightarrow u(z,t)",
            size=36,
            color=INK,
            width=11.2,
        ).next_to(steps, DOWN, buff=0.72)
        takeaways = VGroup(
            txt("一致部分：动态系数、孔弹波模、六阶边界系统、Sommerfeld 结构。", size=24, color=INK, width=10.5),
            txt("主动取舍：Ricker 压力源、零偏移、单界面、Pf 共震项默认关闭。", size=24, color=RED, width=10.5),
            txt("写作重点：始终说明哪些是 Schakel 理论，哪些是 RT-SE 项目映射。", size=24, color=INK, width=10.5),
        ).arrange(DOWN, buff=0.18).next_to(key_formula, DOWN, buff=0.55)
        self.play(FadeIn(header), FadeIn(steps, lag_ratio=0.1), run_time=1.0)
        self.play(ShowCreation(arrows, lag_ratio=0.12), FadeIn(key_formula), run_time=1.2)
        self.play(FadeIn(takeaways, lag_ratio=0.15), run_time=1.0)
        self.active = VGroup(header, steps, arrows, key_formula, takeaways)
        self.hold("summary")
