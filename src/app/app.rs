use log::{info, error};
use clap::{Parser, Subcommand};
use super::cmd::{kh, ckh, ss, ss_batch};
use crate::utils::*;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
pub struct CliArgs {
    #[command(subcommand)]
    pub command: Cmd
}

#[derive(Subcommand, Debug)]
pub enum Cmd {
    Kh(kh::Args),
    Ckh(ckh::Args),
    SS(ss::Args),
    SSBatch(ss_batch::Args),
}

impl Cmd { 
    fn debug(&self) -> bool { 
        match self { 
            Cmd::Kh(args)  => args.debug,
            Cmd::Ckh(args) => args.debug,
            Cmd::SS(args)  => args.debug,
            Cmd::SSBatch(args)  => args.debug
        }
    }
}

pub struct App {}

impl App { 
    pub fn new() -> Self { 
        App {}
    }

    pub fn run(&self) -> Result<String, i32> { 
        let args = CliArgs::parse();

        if args.command.debug() { 
            self.init_logger();
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

    fn init_logger(&self) {
        use simplelog::*;
        TermLogger::init(
            LevelFilter::Info,
            Config::default(),
            TerminalMode::Mixed,
            ColorChoice::Auto
        ).unwrap()
    }

    fn dispatch(&self, args: &CliArgs) -> Result<String, Box<dyn std::error::Error>> { 
        guard_panic(||
            match &args.command { 
                Cmd::Kh(args)      => kh::run(args),
                Cmd::Ckh(args)     => ckh::run(args),
                Cmd::SS(args)      => ss::run(args),
                Cmd::SSBatch(args) => ss_batch::run(args),
            }
        )
    }
}