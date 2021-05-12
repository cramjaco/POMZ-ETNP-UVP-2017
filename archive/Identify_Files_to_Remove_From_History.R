# These are files of unused data that I leave to Andrew McDonnell to make public. Please contact the authors if you want any of them.

source("UVP_2017_library.R")

uvpMeta <- read_tsv("data/uvpdata/export_detailed_20190304_23_14_Export_metadata_summary.tsv") %>% rename(time = `yyyy-mm-dd hh:mm`) %>% arrange(time)

# We don't want particle files form stations other than 016
uvpMeta %>% filter(Site != "016") %>% pull(`Particle filename`) %>% paste0("data/uvpdata/", .) %>% write("FilesToRemoveFromHistory.txt")

# and we don't want zooplankton files
uvpMeta %>% pull(`Plankton filename`) %>% paste0("data/uvpdata/", .) %>% write("FilesToRemoveFromHistory.txt", append = TRUE)

## In a cloned repo, we will run, in the terminal
# https://htmlpreview.github.io/?https://github.com/newren/git-filter-repo/blob/docs/html/git-filter-repo.html
# git filter-repo --invert-paths --paths-from-file /tmp/files-i-dont-want-anymore.txt
# git filter-repo --invert-paths --paths-from-file FilesToRemoveFromHistory.txt