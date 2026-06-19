from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "tmp" / "schakel2011_sommerfeld_explainer"
ASSETS = ROOT / "hyperframes" / "schakel2011_sommerfeld_explainer" / "assets"
VOICE = "zh-CN-XiaoxiaoNeural"
RATE = "+0%"
EDGE_TTS = "edge-tts.cmd" if os.name == "nt" else "edge-tts"
FFMPEG = "ffmpeg.exe" if os.name == "nt" else "ffmpeg"
FFPROBE = "ffprobe.exe" if os.name == "nt" else "ffprobe"


SEGMENTS = [
    {
        "key": "title",
        "pause": 1.0,
        "sections": "title, 1",
        "title": "总览：这段代码在解决什么",
        "text": "这条动画讲的是 schakel 二零一一年 Sommerfeld 正演代码。它的任务不是重新发明震电理论，而是把反应输运给出的孔隙结构和化学变化，接到 Schakel 系列论文的界面电磁响应模型上，最后合成接收点随深度变化的波形。",
    },
    {
        "key": "two_times",
        "pause": 1.0,
        "sections": "1, 2, 10",
        "title": "两条时间轴",
        "text": "先抓住一个很容易混淆的点。溶蚀时间是材料慢慢演化的时间，可以是秒到小时。波形时间是声源激发后的传播时间，通常只有微秒。代码用 Time_s 选择每一个材料状态，再在这个固定状态下计算一次频域响应和时间域波形。相位约定采用时间项 exp i omega t，空间项 exp 负 i k 点 x；界面上方是流体侧 z 小于零，下方是多孔介质侧 z 大于零，复垂向波数分支要保证沿传播方向衰减。",
    },
    {
        "key": "rt_inputs",
        "pause": 1.0,
        "sections": "3",
        "title": "反应输运参数入口",
        "text": "每一个溶蚀时刻，代码读取孔隙率、渗透率、曲折度和氢离子浓度。渗透率从毫达西换成平方米，氢离子浓度换成摩尔每升并计算 pH。电解质浓度和 zeta 电位是项目级映射假设，需要和论文里的理论公式分开表述。",
    },
    {
        "key": "dynamic_bridge",
        "pause": 1.2,
        "sections": "4",
        "title": "动态电动桥梁",
        "text": "接下来进入 Schakel 和 Smeulders 二零一零附录 A。动态渗透率描述孔隙流体在不同频率下是否还能跟上压力振荡。Debye 长度和 zeta 电位进入电动耦合系数。动态电导率和有效介电常数则告诉电磁场在多孔介质中如何响应。这里要特别区分，多孔介质内部的动态体电导率 sigma omega，不等于上覆自由流体电导率 sigma fl，后者只进入流体侧电磁慢度。",
    },
    {
        "key": "elastic_modes",
        "pause": 1.0,
        "sections": "5",
        "title": "孔弹波模",
        "text": "材料参数还要进入 Biot 孔弹框架。代码先计算 A、Q、R 和 P，再计算频率相关密度以及电动耦合修正项。最后得到快 P 波、慢 P 波、TM 波和 SV 波四类波模，以及把机械势和电磁势联系起来的 beta 与 alpha 幅值比。",
    },
    {
        "key": "interface_matrix",
        "pause": 1.2,
        "sections": "6",
        "title": "六个边界条件",
        "text": "界面处真正发生转换。流体侧有入射声波、反射声波和反射电磁波；多孔介质侧有快 P、慢 P、TM 和 SV 波。open pore 条件要求法向体积位移连续、压力连续、骨架应力满足界面条件，同时切向磁场和电场连续。六个条件组成一个六阶线性系统。",
    },
    {
        "key": "coefficients",
        "pause": 1.0,
        "sections": "6, 7",
        "title": "从 2010 系数到 2011 系数",
        "text": "六阶系统解出的核心量是流体侧反射电磁系数、多孔侧 TM 透射系数和快 P 透射系数。但二零一零论文的系数是位移归一化，二零一一年 Sommerfeld 正演需要压力归一化。所以代码把反射电磁系数除以流体密度和角频率平方，多孔侧电势项还要乘上对应的 alpha 幅值比。",
    },
    {
        "key": "source",
        "pause": 1.0,
        "sections": "8.1, 8.2",
        "title": "声源谱和圆形活塞方向性",
        "text": "二零一一年论文用实验压力记录作为声源谱。当前代码默认用因果 Ricker 压力波形，并加平滑带通窗和起始 ramp。圆形换能器的方向性由一阶 Bessel 函数描述，放进 Sommerfeld 积分后，会和零阶 Bessel 函数一起控制不同角度的贡献。",
    },
    {
        "key": "sommerfeld",
        "pause": 1.3,
        "sections": "8.3",
        "title": "Sommerfeld 积分",
        "text": "流体侧反射电势不是一个单一角度的平面波，而是把所有径向波数贡献叠加起来。代码把积分路径拆成两段：实角分支代表传播波，倏逝分支代表沿界面快速衰减但近场重要的贡献。反射电磁系数随频率和角度变化，所以积分自然给出接收深度相关的界面电磁波形。",
    },
    {
        "key": "porous_side",
        "pause": 1.0,
        "sections": "9",
        "title": "多孔侧为什么是简化模型",
        "text": "Schakel 二零一一年 JAP 的完整实验模型包含有限厚度样品、后界面反射和多次响应。当前代码是单界面模型，默认只保留前界面的 TM interface EM 项。快 P 共震电势项作为诊断开关存在，但默认关闭，用来突出界面响应而不是样品多次反射。",
    },
    {
        "key": "time_domain",
        "pause": 1.0,
        "sections": "10, 11, 12",
        "title": "回到时间域与诊断",
        "text": "频域响应算完以后，代码只显式积分正频率，再取两倍实部恢复时间波形。界面到时 T 零由声源到界面的距离和流体声速决定。接收点分为流体侧和多孔侧，z 等于零的界面行只是绘图插值，不参与峰值和诊断统计。诊断包括积分收敛、界面上下极性反转，以及 T 零前后能量比。",
    },
    {
        "key": "summary",
        "pause": 1.5,
        "sections": "13, 14, 15",
        "title": "论文写作中的一句话机制链",
        "text": "最后把公式压缩成一句机制链：溶蚀改变孔隙连通性和流体化学状态；这些变化进入动态渗透率、电动耦合系数和动态电导率；界面边界条件把材料变化转成压力归一化转换系数；Sommerfeld 积分再把角度和频率相关的转换合成为接收点波形。这就是反应输运到 interface EM response 的证据链。",
    },
]


def run(cmd: list[str]) -> None:
    subprocess.run(cmd, cwd=ROOT, check=True)


def ffprobe_duration(path: Path) -> float:
    result = subprocess.run(
        [
            FFPROBE,
            "-v",
            "error",
            "-show_entries",
            "format=duration",
            "-of",
            "default=noprint_wrappers=1:nokey=1",
            str(path),
        ],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return float(result.stdout.strip())


def main() -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    ASSETS.mkdir(parents=True, exist_ok=True)

    concat_lines: list[str] = []
    timeline = []
    current = 0.0

    for index, segment in enumerate(SEGMENTS, start=1):
        key = segment["key"]
        text_path = WORK / f"{index:02d}_{key}.txt"
        audio_path = ASSETS / f"{index:02d}_{key}.mp3"
        silence_path = ASSETS / f"{index:02d}_{key}_pause.mp3"
        text_path.write_text(segment["text"], encoding="utf-8")

        run(
        [
                EDGE_TTS,
                "--voice",
                VOICE,
                "--rate",
                RATE,
                "--file",
                str(text_path),
                "--write-media",
                str(audio_path),
            ]
        )

        pause = float(segment["pause"])
        run(
            [
                FFMPEG,
                "-y",
                "-f",
                "lavfi",
                "-i",
                "anullsrc=r=24000:cl=mono",
                "-t",
                f"{pause:.3f}",
                "-q:a",
                "9",
                "-acodec",
                "libmp3lame",
                str(silence_path),
            ]
        )

        voice_duration = ffprobe_duration(audio_path)
        end = current + voice_duration + pause
        timeline.append(
            {
                "key": key,
                "title": segment["title"],
                "sections": segment["sections"],
                "start": round(current, 3),
                "voice_end": round(current + voice_duration, 3),
                "end": round(end, 3),
                "duration": round(voice_duration + pause, 3),
                "text": segment["text"],
            }
        )
        current = end
        concat_lines.append(f"file '{audio_path.as_posix()}'")
        concat_lines.append(f"file '{silence_path.as_posix()}'")

    concat_path = WORK / "concat_audio.txt"
    concat_path.write_text("\n".join(concat_lines) + "\n", encoding="utf-8")
    run(
        [
            FFMPEG,
            "-y",
            "-f",
            "concat",
            "-safe",
            "0",
            "-i",
            str(concat_path),
            "-ar",
            "24000",
            "-ac",
            "1",
            "-c:a",
            "libmp3lame",
            "-q:a",
            "4",
            str(ASSETS / "voiceover.mp3"),
        ]
    )

    (ASSETS / "timeline.json").write_text(
        json.dumps(
            {
                "voice": VOICE,
                "rate": RATE,
                "total_duration": round(current, 3),
                "segments": timeline,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    (WORK / "coverage_table.json").write_text(
        json.dumps(
            [
                {
                    "key": s["key"],
                    "source_sections": s["sections"],
                    "title": s["title"],
                }
                for s in SEGMENTS
            ],
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    print("voiceover=hyperframes/schakel2011_sommerfeld_explainer/assets/voiceover.mp3")
    print("timeline=hyperframes/schakel2011_sommerfeld_explainer/assets/timeline.json")
    print(f"duration={current:.3f}")


if __name__ == "__main__":
    main()
