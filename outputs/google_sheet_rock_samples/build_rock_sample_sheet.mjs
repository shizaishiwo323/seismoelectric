import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const outDir = "C:/Users/imgw/Documents/Codex/震电模拟/outputs/google_sheet_rock_samples";
const photoDir = "C:/Users/imgw/Downloads/碳酸岩样品";

const samples = [
  {
    id: "M14/24",
    k: 124.3,
    phi: 38.2,
    frontSize: "待量取",
    sideSize: "待量取",
    frontPhoto: "24-正.jpg",
    sidePhoto: "24-侧.jpg",
  },
  {
    id: "M14/25",
    k: 107.6,
    phi: 37.4,
    frontSize: "待量取",
    sideSize: "待量取",
    frontPhoto: "25-正.jpg",
    sidePhoto: "25-侧.jpg",
  },
  {
    id: "M14/26",
    k: 59.7,
    phi: 34.2,
    frontSize: "待量取",
    sideSize: "待量取",
    frontPhoto: "26-正.jpg",
    sidePhoto: "26-侧.jpg",
  },
  {
    id: "M14/34",
    k: 108.7,
    phi: 35.6,
    frontSize: "待量取",
    sideSize: "待量取",
    frontPhoto: "34-正.jpg",
    sidePhoto: "34-侧.jpg",
  },
];

async function imageDataUrl(filename) {
  const filePath = path.join(photoDir, filename);
  const bytes = await fs.readFile(filePath);
  return `data:image/jpeg;base64,${bytes.toString("base64")}`;
}

await fs.mkdir(outDir, { recursive: true });

const workbook = Workbook.create();
const sheet = workbook.worksheets.add("碳酸岩样品");
sheet.showGridLines = false;

sheet.getRange("A1:G1").values = [[
  "样品编号",
  "渗透率（沿钻孔方向）",
  "孔隙率",
  "正面尺寸",
  "侧面尺寸",
  "正面照片",
  "侧面照片",
]];

sheet.getRange(`A2:E${samples.length + 1}`).values = samples.map((sample) => [
  sample.id,
  sample.k,
  sample.phi / 100,
  sample.frontSize,
  sample.sideSize,
]);

sheet.getRange("A1:G1").format = {
  fill: "#1F4E79",
  font: { bold: true, color: "#FFFFFF" },
  horizontalAlignment: "center",
  verticalAlignment: "middle",
  wrapText: true,
};
sheet.getRange(`A2:E${samples.length + 1}`).format = {
  verticalAlignment: "middle",
  wrapText: true,
};
sheet.getRange(`A2:A${samples.length + 1}`).format.horizontalAlignment = "center";
sheet.getRange(`B2:C${samples.length + 1}`).format.horizontalAlignment = "right";
sheet.getRange(`D2:E${samples.length + 1}`).format.horizontalAlignment = "center";
sheet.getRange(`A1:G${samples.length + 1}`).format.borders = {
  preset: "all",
  style: "thin",
  color: "#D9E2F3",
};
sheet.getRange("B2:B5").format.numberFormat = "0.0";
sheet.getRange("C2:C5").format.numberFormat = "0.0%";

sheet.getRange("A:A").format.columnWidthPx = 95;
sheet.getRange("B:B").format.columnWidthPx = 155;
sheet.getRange("C:C").format.columnWidthPx = 85;
sheet.getRange("D:E").format.columnWidthPx = 95;
sheet.getRange("F:G").format.columnWidthPx = 190;
sheet.getRange("1:1").format.rowHeightPx = 48;
sheet.getRange(`2:${samples.length + 1}`).format.rowHeightPx = 142;
sheet.freezePanes.freezeRows(1);

for (let i = 0; i < samples.length; i += 1) {
  const row = i + 1;
  const sample = samples[i];
  sheet.images.add({
    dataUrl: await imageDataUrl(sample.frontPhoto),
    anchor: { from: { row, col: 5, rowOffsetPx: 8, colOffsetPx: 10 }, extent: { widthPx: 170, heightPx: 124 } },
  });
  sheet.images.add({
    dataUrl: await imageDataUrl(sample.sidePhoto),
    anchor: { from: { row, col: 6, rowOffsetPx: 8, colOffsetPx: 10 }, extent: { widthPx: 170, heightPx: 124 } },
  });
}

const notes = workbook.worksheets.add("来源说明");
notes.showGridLines = false;
notes.getRange("A1:B6").values = [
  ["项目", "说明"],
  ["孔隙率/渗透率来源", "C:/Users/imgw/Downloads/碳酸岩样品/孔隙率.jpg"],
  ["照片来源", "C:/Users/imgw/Downloads/碳酸岩样品"],
  ["渗透率单位", "mD，按孔隙率.jpg 中 K(mD) 手写列录入"],
  ["孔隙率单位", "%，表格内以百分比格式显示"],
  ["尺寸列", "源照片含尺规但未给出明确数值，已保留为待量取，避免误填"],
];
notes.getRange("A1:B1").format = {
  fill: "#1F4E79",
  font: { bold: true, color: "#FFFFFF" },
};
notes.getRange("A:B").format.columnWidthPx = 260;
notes.getRange("A1:B6").format.borders = { preset: "all", style: "thin", color: "#D9E2F3" };

const preview = await workbook.render({
  sheetName: "碳酸岩样品",
  range: "A1:G5",
  scale: 1,
  format: "png",
});
await fs.writeFile(path.join(outDir, "rock_sample_sheet_preview.png"), new Uint8Array(await preview.arrayBuffer()));

const errorScan = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 100 },
  summary: "final formula error scan",
});
console.log(errorScan.ndjson);

const xlsx = await SpreadsheetFile.exportXlsx(workbook);
await xlsx.save(path.join(outDir, "碳酸岩样品汇总.xlsx"));
console.log(path.join(outDir, "碳酸岩样品汇总.xlsx"));
