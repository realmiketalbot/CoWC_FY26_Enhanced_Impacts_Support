// AGU 2026 poster wireframe: 42 in wide x 36 in tall, single slide.
// Figures are placed at exactly the size they were exported at (300 dpi), so
// their text prints at the AGU 2025 sizes. The decision tree is drawn with
// native, editable shapes.
const pptxgen = require("pptxgenjs");
const path = require("path");

const FIG = "/Users/mtalbot/Desktop/GitHub/CoWC_FY26_Enhanced_Impacts_Support/_documentation/poster_figures";
const OUT = process.argv[2];

const GREEN = "1E4D2B", GOLD = "C8C372", GREY = "F2F2F2", DARK = "222222", MUTED = "595959";
const C_RAIN = "B2182B", C_SNOW = "2166AC", C_ROS = "9970AB", C_UNC = "CCCCCC", C_DIAM = "FFF5C2";
const FONT = "Arial";

const pres = new pptxgen();
pres.defineLayout({ name: "POSTER_42x36", width: 42, height: 36 });
pres.layout = "POSTER_42x36";
pres.title = "AGU 2026 poster wireframe";
const s = pres.addSlide();
s.background = { color: "FFFFFF" };

// ---- geometry (must match PANEL_W / WIDE_W in poster_figures.R) ----
const M = 0.6, G = 0.5;
const W = 9.825, W2 = 20.15;                       // one column, two columns
const colX = (i) => M + i * (W + G);
const TOP = 5.0, BOT = 36 - M, HDR = 1.0, HGAP = 0.15;

const text = (t, x, y, w, h, o = {}) =>
  s.addText(t, Object.assign({ x, y, w, h, fontFace: FONT, fontSize: 24, color: DARK, valign: "top",
                               margin: 0.1, isTextBox: true }, o));
const header = (t, x, y, w) => {
  s.addShape(pres.shapes.RECTANGLE, { x, y, w, h: HDR, fill: { color: GREEN }, line: { color: GREEN } });
  text(t, x, y, w, HDR, { fontSize: 30, bold: true, color: "FFFFFF", valign: "middle", margin: 0.25 });
  return y + HDR + HGAP;
};
const placeholder = (label, x, y, w, h) => {
  s.addShape(pres.shapes.RECTANGLE, { x, y, w, h, fill: { color: GREY },
                                      line: { color: "7F7F7F", width: 2, dashType: "dash" } });
  text(label, x, y, w, h, { fontSize: 22, italic: true, color: MUTED, align: "center", valign: "middle" });
};
// Place a figure at its exported size (w x h inches)
const fig = (file, x, y, w, h) => { s.addImage({ path: path.join(FIG, file), x, y, w, h }); return y + h; };
// Lists are inset from the section boxes by LIST_PAD so bullets don't sit flush with the box edge.
const LIST_PAD = 0.5;
const bullets = (items, x, y, w, h, size = 24) =>
  text(items.map((t, i) => ({ text: t, options: { bullet: true, breakLine: i < items.length - 1, paraSpaceAfter: 6 } })),
       x + LIST_PAD, y, w - LIST_PAD, h, { fontSize: size });

// ---- header band ----
placeholder("AGU26 logo", M, 0.6, 4.5, 3.4);
placeholder("CSU logo", 42 - M - 4.5, 0.6, 4.5, 3.4);
text("A Regional Mixed-Population Flood Frequency Analysis: Insights from Colorado",
     M + 5, 0.35, 42 - 2 * M - 10, 2.3, { fontSize: 58, bold: true, color: GREEN, align: "center", valign: "middle" });
text("Michael Talbot and Frances V. Davenport", M + 5, 2.65, 42 - 2 * M - 10, 0.9,
     { fontSize: 40, align: "center", valign: "middle" });
text("Civil and Environmental Engineering, Colorado State University  |  Preliminary results",
     M + 5, 3.45, 42 - 2 * M - 10, 0.75, { fontSize: 28, color: MUTED, align: "center", valign: "middle" });
s.addShape(pres.shapes.LINE, { x: M, y: 4.5, w: 42 - 2 * M, h: 0, line: { color: GOLD, width: 6 } });

// ---- column 1: motivation + methods ----
let x = colX(0), y = TOP;
y = header("Colorado floods have more than one cause", x, y, W);
bullets([
  "Flood frequency analysis (FFA) assumes peaks come from one, identically distributed population",
  "Colorado floods arise from rainfall, snowmelt and rain-on-snow, which may follow different distributions",
  "Q: Does a mixed-population FFA [1] change design-flood estimates, at which return periods, and where?",
], x, y, W, 3.6, 22);
y += 3.6;
y = fig("02b_peak_timing_by_type.png", x, y, W, 1200 / 300) + G;

y = header("Each peak classified by rain fraction", x, y, W);
bullets([
  "215 USGS gages [2], WY 1986–2023; ≥ 20 years with ≥ 90% May–Sep coverage",
  "Independent peaks, declustered [3, 4]",
  "WLDAS [5] rainfall + snowmelt over a basin-scaled antecedent window [6]",
  "Poisson–GPD [7] per flood type (L-moments [8, 9]); mixture = product of annual-maximum CDFs, vs one fit to all peaks",
], x, y, W, 4.6, 20);
y += 4.6;

// Native decision tree (editable). Current Script 6b rules.
const tree = (x0, y0, w) => {
  const node = (shape, label, cx, cy, nw, nh, fill, fc = DARK, bold = false) => {
    s.addShape(shape, { x: cx - nw / 2, y: cy - nh / 2, w: nw, h: nh, fill: { color: fill },
                        line: { color: "404040", width: 1.5 } });
    text(label, cx - nw / 2, cy - nh / 2, nw, nh, { fontSize: 16, color: fc, bold, align: "center", valign: "middle", margin: 0.05 });
  };
  const arrow = (x1, y1, x2, y2, lab, lx, ly) => {
    s.addShape(pres.shapes.LINE, { x: Math.min(x1, x2), y: Math.min(y1, y2), w: Math.abs(x2 - x1) || 0.001,
                                   h: Math.abs(y2 - y1) || 0.001, flipH: x2 < x1,
                                   line: { color: "404040", width: 1.5, endArrowType: "triangle" } });
    if (lab) text(lab, lx, ly, 0.7, 0.4, { fontSize: 14, color: MUTED, margin: 0 });
  };
  const cx = x0 + w * 0.36, ox = x0 + w * 0.82, dw = 4.2, dh = 1.25, ew = 2.6, eh = 0.75;
  const r1 = y0 + 0.6, r2 = r1 + 1.75, r3 = r2 + 1.75, r4 = r3 + 1.45;
  node(pres.shapes.DIAMOND, "Enough WLDAS data and\nflux ≥ 0.1 mm?", cx, r1, dw, dh, C_DIAM);
  node(pres.shapes.OVAL, "Unclassified", ox, r1, ew, eh, C_UNC);
  arrow(cx + dw / 2, r1, ox - ew / 2, r1, "No", cx + dw / 2 + 0.1, r1 - 0.42);
  node(pres.shapes.DIAMOND, "Rain fraction ≥ 0.99?", cx, r2, dw, dh, C_DIAM);
  arrow(cx, r1 + dh / 2, cx, r2 - dh / 2, "Yes", cx + 0.1, r1 + dh / 2);
  node(pres.shapes.OVAL, "Rainfall", ox, r2, ew, eh, C_RAIN, "FFFFFF", true);
  arrow(cx + dw / 2, r2, ox - ew / 2, r2, "Yes", cx + dw / 2 + 0.1, r2 - 0.42);
  node(pres.shapes.DIAMOND, "Rain fraction ≥ network\nmedian of mixed (0.31)?", cx, r3, dw, dh, C_DIAM);
  arrow(cx, r2 + dh / 2, cx, r3 - dh / 2, "No", cx + 0.1, r2 + dh / 2);
  node(pres.shapes.OVAL, "Rain-on-Snow", ox, r3, ew, eh, C_ROS, "FFFFFF", true);
  arrow(cx + dw / 2, r3, ox - ew / 2, r3, "Yes", cx + dw / 2 + 0.1, r3 - 0.42);
  node(pres.shapes.OVAL, "Snowmelt", cx, r4, ew, eh, C_SNOW, "FFFFFF", true);
  arrow(cx, r3 + dh / 2, cx, r4 - eh / 2, "No", cx + 0.1, r3 + dh / 2);
  text("2-class: Rain-on-Snow + Snowmelt = Snow-Influenced", x0, r4 + eh / 2 + 0.05, w, 0.45,
       { fontSize: 16, italic: true, color: MUTED, align: "center" });
  return r4 + eh / 2 + 0.55;
};
y = tree(x, y, W) + 0.1;
y = fig("02_composition_by_method.png", x, y, W, 2100 / 300);
console.log("column 1 ends at", y.toFixed(2), "of", BOT);

// ---- columns 2–3: example site + headline result ----
x = colX(1); y = TOP;
y = header("At Tomichi Creek, the mixture rises above the single sample for rare floods", x, y, W2);
y = fig("showcase/09115500_A_fixed_vs_adaptive3.png", x, y, W2, 3000 / 300);
bullets([
  "The depth threshold method [10] misses many peaks, so its mixture falls below observed common floods; the rain fraction method tracks the observations",
  "100-yr mixture / single sample: 1.30 (rain fraction, 3-class) vs 1.16 (depth threshold)",
], x, y + 0.05, W2, 1.4, 22);
y += 1.45 + G;

y = header("Estimates agree for common floods and diverge for rare ones", x, y, W2);
fig("04_ratio_by_return_period.png", x, y, W, 2100 / 300);
y = fig("06_ratio_map_q100.png", x + W + G, y, W, 2100 / 300) + 0.2;
const stats = [["85%", "of sites: mixture > single sample at 100 yr"],
               ["+10%", "median difference at 100 yr (rain fraction, 3-class)"],
               ["~0%", "median difference at 2 yr"]];
const sw = (W2 - 2 * G) / 3, sh = 2.2;
stats.forEach(([big, small], i) => {
  const sx = x + i * (sw + G);
  s.addShape(pres.shapes.RECTANGLE, { x: sx, y, w: sw, h: sh, fill: { color: "EEF2EC" }, line: { color: "EEF2EC" } });
  text(big, sx, y + 0.05, sw, 1.15, { fontSize: 56, bold: true, color: GREEN, align: "center", valign: "middle" });
  text(small, sx, y + 1.2, sw, 0.95, { fontSize: 20, color: DARK, align: "center" });
});
y += sh + 0.15;
text([{ text: "Robust to analysis choices: ", options: { bold: true } },
      { text: "median 100-yr ratio 1.05–1.28 across fitting thresholds and 1.10–1.16 across rain-fraction cutoffs. " },
      { text: "Depends on classification: ", options: { bold: true } },
      { text: "depth threshold 1.16, rain fraction 3-class 1.10, rain fraction 2-class 1.01, so separating rain-on-snow from snowmelt drives most of it." }],
     x, y, W2, BOT - y, { fontSize: 20 });

// ---- column 4: elevation, flood type, takeaways, next steps ----
x = colX(3); y = TOP;
y = header("Differences peak in high-elevation basins", x, y, W);
y = fig("08b_ratio_by_elevation_band.png", x, y, W, 1800 / 300) + G;
y = header("Rain-on-snow drives more rare floods", x, y, W);
y = fig("09_dominant_type_by_rp.png", x, y, W, 1800 / 300) + G;

y = header("Takeaways", x, y, W);
text([
  { text: "1  Floods at nearly all Colorado gages come from more than one process (91% of sites)", options: { breakLine: true, paraSpaceAfter: 8 } },
  { text: "2  Mixture and single-sample estimates agree for common floods and diverge for rare ones, most strongly above 3,000 m", options: { breakLine: true, paraSpaceAfter: 8 } },
  { text: "3  Separating rain-on-snow from snowmelt matters most. Is that split physically distinct?" },
], x + LIST_PAD, y, W - LIST_PAD, 4.2, { fontSize: 21 });
y += 4.2 + 0.2;
y = header("Next steps", x, y, W);
const qr = 2.4;
bullets([
  "EQD threshold selection [11]; double-bootstrap confidence intervals",
  "Rain-on-snow vs snowmelt diagnostic; goodness-of-fit",
  "Trends: no detectable change in annual-maximum flood type (1986–2023)",
], x, y, W - qr - 0.2, qr, 17);
placeholder("New QR\n(PIN)", x + W - qr, y, qr, qr);
y += qr + G;

// References: entries from the Preliminary Proposal's bibliography, in the
// numbered style of the AGU 2025 poster (10.5 pt).
y = header("References", x, y, W);
const REFS = [
  "Yu G, Wright DB, Davenport FV. Diverse physical processes drive upper-tail flood quantiles in the US Mountain West. Geophys Res Lett. 2022;49:e2022GL098855.",
  "Falcone JA. GAGES-II: Geospatial Attributes of Gages for Evaluating Streamflow. U.S. Geological Survey; 2011. doi:10.3133/70046617",
  "USWRC. Guidelines for determining flood flow frequency: Bulletin 17. U.S. Water Resources Council; 1975.",
  "Lang M, Ouarda TBMJ, Bobée B. Towards operational guidelines for over-threshold modeling. J Hydrol. 1999;225:103–117.",
  "Erlingis JM et al. A high-resolution land data assimilation system optimized for the western United States. JAWRA. 2021;57:692–710.",
  "Livneh B et al. Can precipitation intermittency predict flooding? Sci Total Environ. 2024;945:173824.",
  "Beguería S. Uncertainties in partial duration series modelling of extremes related to the choice of the threshold value. J Hydrol. 2005;303:215–230.",
  "Hosking JRM, Wallis JR. Parameter and quantile estimation for the generalized Pareto distribution. Technometrics. 1987;29:339–349.",
  "Asquith W. lmomco: L-moments, censored L-moments, trimmed L-moments, L-comoments, and many distributions. R package; 2026.",
  "Davenport FV, Herrera-Estrada JE, Burke M, Diffenbaugh NS. Flood size increases nonlinearly across the western United States in response to lower snow-precipitation ratios. Water Resour Res. 2020;56:e2019WR025571.",
  "Murphy C, Tawn JA, Varty Z. Automated threshold selection and associated inference uncertainty for univariate extremes. Technometrics. 2025;67:215–224.",
];
text(REFS.map((r, i) => ({ text: `${i + 1}. ${r}`, options: { breakLine: i < REFS.length - 1, paraSpaceAfter: 2 } })),
     x + LIST_PAD, y, W - LIST_PAD, BOT - y - 0.45, { fontSize: 10.5 });
text("Funding: FY26 Colorado Water Center Enhanced Impacts Initiative", x, BOT - 0.4, W, 0.4,
     { fontSize: 14, color: MUTED, valign: "bottom" });

console.log("column 4 ends at", y.toFixed(2), "of", BOT);
pres.writeFile({ fileName: OUT }).then(() => console.log("wrote " + OUT));
