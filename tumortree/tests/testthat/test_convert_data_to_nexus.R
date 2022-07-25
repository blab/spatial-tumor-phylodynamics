context("Tests for convert_sim_to_nexus.R")

library(dplyr)
library(tumortree)

#from home dir when writing tests
#test_cells <- read.csv("tumortree/tests/testdata/test_cells.csv")

test_cells <- read.csv("../testdata/test_cells.csv")

normalized_test_cells <- normalize_locs(test_cells)

alive_test_cells <- test_cells %>%
    dplyr::filter(deathdate == max(test_cells$deathdate))

sampled_test_cells <- alive_cells

#tests for extract_mutations
sampled_muts <- sort(unique(as.numeric(unlist(purrr::map(1:nrow(sampled_cells),
                                             function(i) extract_mutations(sampled_cells$mutations[i]))))))
max_mut_number <- max(as.numeric(unlist(purrr::map(1:nrow(test_cells),
                                                   function(i) extract_mutations(test_cells$mutations[i])))),
                      na.rm = TRUE)

row_one_muts <- extract_mutations(test_cells$sequence[1])

test_that('Extracted sequence only contains nucleotides', {
    expect_true(all(row_one_muts %in% c("A", "C", "G", "T")))
})


test_that('Length of extracted sequence is equal to the number of mutations in cells', {
    expect_equal(length(row_one_muts), max_mut_number)
})

# sim_cells <- test_cells %>%
#     normalize_locs %>%
#     dplyr::arrange(index)
#
# endpoint <- max(sim_cells$deathdate)
# sampled_cells <- sim_cells %>%
#     dplyr::filter(deathdate == endpoint)
#
# sampled_muts <- sampled_cells$mutations
#
# all_sampled_muts <- sort(as.integer(unique(unlist(purrr::map(sampled_muts, function(mut_row) extract_mutations(mut_row))))))
# ## extract sequence data
#
#
# sampled_cells$sequence_collapsed <- purrr::map_chr(sampled_cells$sequence, function(seq) extract_mutations(seq, sampled_muts = all_sampled_muts))
#

