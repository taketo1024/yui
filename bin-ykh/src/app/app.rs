use log::info;
use clap::{Parser, Subcommand};

use super::cmd::{ckh, ckhi, kh, khi};
use super::utils::*;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
pub struct CliArgs {
    #[command(subcommand)]
    pub command: Cmd
}

#[derive(Subcommand, Debug)]
#[clap(rename_all="lower")]
pub enum Cmd {
    CKh(ckh::Args),
    Kh(kh::Args),
    CKhI(ckhi::Args),    
    KhI(khi::Args),    
}

impl CliArgs { 
    fn log_level(&self) -> log::LevelFilter { 
        use log::LevelFilter::*;
        let level = match &self.command { 
            Cmd::CKh(args)  => args.log,
            Cmd::Kh(args)   => args.log,
            Cmd::CKhI(args) => args.log,
            Cmd::KhI(args)  => args.log,
        };
        match level {
            1 => Info,
            2 => Debug,
            3 => Trace,
            _ => Off,
        }
    }
}

pub struct App {
    pub args: CliArgs
}

impl App { 
    pub fn new() -> Self { 
        let args = CliArgs::parse();
        App { args }
    }

    pub fn run(&self) -> Result<String, Box<dyn std::error::Error>> { 
        self.init_logger();

        info!("args: {:?}", self.args);
        info!("int-type: {}", std::any::type_name::<super::utils::dispatch::Int>());

        let (res, time) = measure(||
            self.dispatch()
        );

        info!("time: {:?}", time);

        res
    }

    fn init_logger(&self) {
        let l = self.args.log_level();
        yui::util::log::init_simple_logger(l).unwrap()
    }

    fn dispatch(&self) -> Result<String, Box<dyn std::error::Error>> { 
        guard_panic(||
            match &self.args.command { 
                Cmd::CKh(args)  => ckh::dispatch(args),
                Cmd::Kh(args)   => kh::dispatch(args),
                Cmd::CKhI(args) => ckhi::dispatch(args),
                Cmd::KhI(args)  => khi::dispatch(args),
            }
        )
    }
}