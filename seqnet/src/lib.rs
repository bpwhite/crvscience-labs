//! `bio-suite-seqnet`: Pure-Rust 3-Stage Hybrid Sequence Analysis Neural Network
//! for multi-locus phylogenetics, multi-database SQLite caching with soft-quarantine
//! priors, species-per-genus chimeric OTU construction, and clean partitioned NEXUS export.

pub mod align_nexus;
pub mod config;
pub mod db;
pub mod entrez;
pub mod gui;
pub mod neural_tensor;
pub mod neurons;
pub mod otu;
pub mod region;
