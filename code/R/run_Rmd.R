



rmarkdown::render("code/Rmd/main_figs.Rmd",
                  knit_root_dir="../..",
                  output_dir="output",
                  output_format="html_document"
)

rmarkdown::render("code/Rmd/supp_figs.Rmd",
                  knit_root_dir="../..",
                  output_dir="output",
                  output_format="word_document"
)
