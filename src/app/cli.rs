use clap::{Parser, Subcommand};
use super::cmd::{kh, ckh, ss};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
pub struct CliArgs {
    #[command(subcommand)]
    pub command: Cmd,

    #[arg(long, default_value_t = false)]
    pub debug: bool
}

#[derive(Subcommand, Debug)]
pub enum Cmd {
    Kh(kh::Args),
    Ckh(ckh::Args),
    SS(ss::Args),
    SSBatch(ss::BatchArgs)
}
