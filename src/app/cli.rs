use clap::{Parser, Subcommand};
use crate::utils::CType;

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
    Kh {
        name: String,

        #[arg(short, long)]
        link: Option<String>,

        #[arg(short, long, default_value = "0")]
        c_value: String,

        #[arg(short = 't', long, default_value = "z")]
        c_type: CType,

        #[arg(short, long)]
        mirror: bool,

        #[arg(short, long)]
        reduced: bool,

        #[arg(short, long)]
        bigraded: bool,
    },

    Ckh {
        name: String,

        #[arg(short, long)]
        link: Option<String>,

        #[arg(short, long, default_value = "0")]
        c_value: String,

        #[arg(short = 't', long, default_value = "z")]
        c_type: CType,

        #[arg(short, long)]
        mirror: bool,

        #[arg(short, long)]
        reduced: bool,

        #[arg(short = 'a', long)]
        with_alpha: bool,
    },

    SS {
        name: String,

        #[arg(short, long)]
        link: Option<String>,

        #[arg(short, long)]
        c_value: String,

        #[arg(short = 't', long, default_value = "z")]
        c_type: CType,
    },

    SSBatch {
        #[arg(short, long)]
        c_value: String,
        
        #[arg(short = 't', long, default_value = "z")]
        c_type: CType,

        #[arg(short, long)]
        data: String,

        #[arg(short, long)]
        output: Option<String>
    },
}
