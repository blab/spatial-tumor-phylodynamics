context("Tests for tumortree")


library(dplyr)
library(phytools)
library(ape)


test_cells <- read.csv("../testdata/test_cells.csv")
normalized_test_cells <- normalize_locs(test_cells)

#tests for normalize_locs()
test_that('Means of normalized coordinates are 0', {
  expect_true(mean(normalized_test_cells$norm_locx == 0) &&
                mean(normalized_test_cells$norm_locy == 0))
})


test_that('Normalized x values are as expected from last version', {
  expect_equal(normalized_test_cells$norm_locx,
               c(-1,0,-1,1,0,0,-1,-2,-1,0,0,-2,-1,1,0,-2,-1 ,2,1))
})

test_that('Normalized y values are as expected from last version', {
  expect_equal(normalized_test_cells$norm_locy,
               c(0,-1 ,0,-1,-1,1,0,-1,0,0,1,0,0,1,0,1,0,0,1))
})

#tests for filter_alive_cells()
alive_test_cells <- filter_alive_cells(normalized_test_cells)

test_that('All cell deathdates are simulation endpoint', {
  expect_true(all(alive_test_cells$deathdate == max(test_cells$deathdate)))
})


test_that('Everything excluded from the alive cells have died', {
  expect_true(all(dplyr::anti_join(normalized_test_cells, alive_test_cells)$deathdate < max(test_cells$deathdate)))

})

#tests for process_and_sample_cells()
sampled_test_cells <- process_and_sample_cells(normalized_test_cells, n = 5)

test_that('Correct number of cells sampled', {
  sampled_test_cells <- process_and_sample_cells(normalized_test_cells, n = 5)
  expect_true(nrow(sampled_test_cells) == 5)
})


sampled_test_cells <- read.csv("../testdata/sampled_test_cells.csv")
#tests for isolate_muts()
sampled_muts <- isolate_muts(mutations_column = sampled_test_cells$mutations)

test_that('Sampled mutations extracted are the same as expected from last version', {
  expect_equivalent(sampled_muts, c(1,2,4,5,6,7,8,9,12,13,16,17,18))
})

#tests for define_mut_presence_absence()

test_that('presence_absence df has expected number of rows', {
  expect_equal(nrow(sampled_test_cells), nrow(define_mut_presence_absence(sampled_test_cells)))
})

test_that('presence_absence df has expected number of columns', {
  expect_equal(length(sampled_muts), ncol(define_mut_presence_absence(sampled_test_cells)))
})

test_that('First and third mutation position is as expected', {
  expect_true(define_mut_presence_absence(sampled_test_cells)[,1] ==  c(1,0,0,0,0) &&
                define_mut_presence_absence(sampled_test_cells)[,3] == c(1,0,0,0,0))
})

#tests for convert_nodes_to_string()
test_that("Binary tree string is formed", {
  nwk_list <- convert_nodes_to_string(sampled_test_cells$index, normalized_test_cells, branch_unit = "time")
  expect_true(ape::is.binary.phylo(phytools::read.newick(text = nwk_list$tree.text)))
})

test_that("All samples are leaves of tree", {
  nwk_list <- convert_nodes_to_string(sampled_test_cells$index, normalized_test_cells, branch_unit = "time")
  expect_true(length(setdiff(as.character(sampled_test_cells$index), phytools::read.newick(text = nwk_list$tree.text)$tip.label)) == 0)
})

test_that("Conversion will not work if given non-existing cell", {
  expect_error(convert_nodes_to_string(c(sampled_test_cells$index,1002), normalized_test_cells, branch_unit = "time"))
})




