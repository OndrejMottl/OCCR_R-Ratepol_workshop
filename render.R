#----------------------------------------------------------#
#
#
#                 OCCD R-Ratepol workshop
#
#                       Render
#
#
#                      O. Mottl
#                        2023
#
#----------------------------------------------------------#
# Render all files to both `.dm` and `.html` files`

library(quarto)
library(rmarkdown)
library(here)

# helper function o render fie in both formats
render_md_and_html <- function(
    file_name,
    file_name_html = file_name,
    toc = FALSE,
    code_download = FALSE) {
  quarto::quarto_render(
    input = here::here(
      paste0(file_name, ".qmd")
    ),
    execute_dir = here::here(),
    output_format = "gfm",
    output_file = paste0(file_name, ".md")
  )

  rmarkdown::render(
    input = here::here(
      paste0(file_name, ".qmd")
    ),
    output_format = rmarkdown::html_document(
      code_download = code_download,
      toc = toc
    ),
    output_file = here::here(
      paste0("docs/", file_name_html, ".html")
    )
  )
}

# README -----
render_md_and_html(file_name = "README", file_name_html = "index")

# workshop_info -----
render_md_and_html(file_name = "workshop_info")

# pre_workshop -----
render_md_and_html(file_name = "pre_workshop")

# step_by_step_guide -----
render_md_and_html(
  file_name = "part_1_geochemical_data",
  code_download = TRUE,
  toc = TRUE
)

render_md_and_html(
  file_name = "part_2_xrf_data",
  code_download = TRUE,
  toc = TRUE
)
