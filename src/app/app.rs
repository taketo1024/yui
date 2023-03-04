use log::{info, error};
use clap::Parser;
use super::cli::{CliArgs, Cmd};
use super::cmd;
use crate::utils::*;

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

        let (res, time) = measure(|| guard_panic(||
            match args.command { 
                Cmd::Kh(args)  => cmd::kh::run(args),
                Cmd::Ckh(args) => cmd::ckh::run(args),
                Cmd::SS(args)  => cmd::ss::run(args),
                Cmd::SSBatch(args) => cmd::ss::run_batch(args)
            }
        ));

        if let Ok(res) = res { 
            info!("time: {:?}", time);
            Ok(res)
        } else {
            if let Err(e) = res { 
                error!("{}", e);
                eprintln!("\x1b[0;31merror\x1b[0m: {e}");
            }
            Err(1)
        }
    }
}