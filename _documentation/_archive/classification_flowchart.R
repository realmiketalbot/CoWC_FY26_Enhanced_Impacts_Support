library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# ============================================================================
# Flood Type Classification Flowchart
# Adaptive window decision tree (Script 6b)
# Thresholds: RAIN_EM_THRESH = 0.95, SNOW_EM_THRESH = 0.10
# ============================================================================

flowchart <- grViz("
digraph classification {

  graph [layout = dot, rankdir = TB, fontname = 'Helvetica', splines = polyline]

  node [fontname = 'Helvetica', fontsize = 14, margin = 0.3]
  edge [fontname = 'Helvetica', fontsize = 12]

  # ---- Input node ----
  input [label = 'POT peak event\n(n-day antecedent window\nvia search_radius)',
         shape = rectangle, style = filled, fillcolor = '#f0f0f0',
         width = 3.5, penwidth = 1.5]

  # ---- Decision nodes ----
  d1 [label = 'Data sufficient?\nn_days_data >= 0.90 x search_radius\nAND total flux >= 1.0 mm',
      shape = diamond, style = filled, fillcolor = '#fef9c3',
      width = 4.5, height = 2.0, penwidth = 1.5, fixedsize = true]

  d2 [label = 'rain_frac >= 0.95?',
      shape = diamond, style = filled, fillcolor = '#fef9c3',
      width = 2.5, height = 1.5, penwidth = 1.5, fixedsize = true]

  d3 [label = 'rain_frac <= 0.10?',
      shape = diamond, style = filled, fillcolor = '#fef9c3',
      width = 2.5, height = 1.5, penwidth = 1.5, fixedsize = true]

  # ---- Terminal nodes ----
  t_unc  [label = 'Unclassified',
          shape = ellipse, style = filled, fillcolor = '#d9d9d9',
          width = 2, height = 0.9, penwidth = 1.5, fixedsize = true]

  t_rain [label = 'Rainfall',
          shape = ellipse, style = filled, fillcolor = '#2166ac',
          fontcolor = 'white', width = 2, height = 0.9, penwidth = 1.5, fixedsize = true]

  t_ros  [label = 'Rain-on-Snow',
          shape = ellipse, style = filled, fillcolor = '#d7191c',
          fontcolor = 'white', width = 2, height = 0.9, penwidth = 1.5, fixedsize = true]

  t_snow [label = 'Snowmelt',
          shape = ellipse, style = filled, fillcolor = '#4dac26',
          fontcolor = 'white', width = 2, height = 0.9, penwidth = 1.5, fixedsize = true]

  # ---- Edges ----
  input -> d1

  d1 -> t_unc  [taillabel = 'No',  labelangle = 45,  labeldistance = 2.5]
  d1 -> d2     [taillabel = 'Yes', labelangle = -45, labeldistance = 2.5]

  d2 -> t_rain [taillabel = 'Yes', labelangle = 45,  labeldistance = 2.5]
  d2 -> d3     [taillabel = 'No',  labelangle = -45, labeldistance = 2.5]

  d3 -> t_snow [taillabel = 'Yes', labelangle = -45, labeldistance = 2.5]
  d3 -> t_ros  [taillabel = 'No',  labelangle = 45,  labeldistance = 2.5]
}
")

# ============================================================================
# Export to PNG
# ============================================================================
out_path <- "_documentation/classification_flowchart_adaptive.png"

svg_str <- export_svg(flowchart)
rsvg_png(charToRaw(svg_str), file = out_path, width = 1200)

message("Saved: ", out_path)