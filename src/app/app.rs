use log::{info, error};
use clap::{Parser, Subcommand};
use super::cmd::{kh, ckh, ss};
use crate::utils::*;

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

pub struct App {
    debug: bool
}

impl App { 
    pub fn new() -> Self { 
        App { debug: false }
    }

    pub fn run(&mut self) -> Result<String, i32> { 
        let args = CliArgs::parse();

        if args.debug { 
            self.debug = true;
            init_logger();
        }

        info!("args: {:?}", args);
        info!("int-type: {}", std::any::type_name::<crate::utils::dispatch::Int>());

        let (res, time) = measure(||
            self.dispatch(&args)
        );

        let res = res.map_err(|e| { 
            error!("{}", e);
            eprintln!("\x1b[0;31merror\x1b[0m: {e}");
            1 // error code
        });

        info!("time: {:?}", time);

        res
    }

    fn dispatch(&self, args: &CliArgs) -> Result<String, Box<dyn std::error::Error>> { 
        guard_panic(||
            match &args.command { 
                Cmd::Kh(args)      => kh::run(args),
                Cmd::Ckh(args)     => ckh::run(args),
                Cmd::SS(args)      => ss::run(args),
                Cmd::SSBatch(args) => ss::run_batch(args)
            }
        )
    }
}