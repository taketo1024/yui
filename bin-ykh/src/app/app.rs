use log::info;
use clap::{Parser, Subcommand};

use super::cmd::{ckh, ckhi, kh, khi, cc, sl2};
use super::args::*;
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
    CC(cc::Args),
    SL2(sl2::Args),
}

impl CliArgs { 
    fn app_args(&self) -> &dyn AppArgs { 
        match &self.command { 
            Cmd::CKh(args)  => args,
            Cmd::Kh(args)   => args,
            Cmd::CKhI(args) => args,
            Cmd::KhI(args)  => args,
            Cmd::CC(args)   => args,
            Cmd::SL2(args)  => args,
        }
    }

    fn log_level(&self) -> log::LevelFilter { 
        self.app_args().log_level()
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

        info!("args:\n{:#?}", self.args);
        info!("int-type: {}", std::any::type_name::<super::utils::dispatch::Int>());

        let (res, time) = measure(||
            self.dispatch()
        );

        info!("time: {:?}", time);

        res
    }

    fn init_logger(&self) {
        let l = self.args.log_level();
        env_logger::Builder::new().filter_level(l).init();
    }

    fn dispatch(&self) -> Result<String, Box<dyn std::error::Error>> { 
        guard_panic(||
            match &self.args.command { 
                Cmd::CKh(args)  => ckh::dispatch(args),
                Cmd::Kh(args)   => kh::dispatch(args),
                Cmd::CKhI(args) => ckhi::dispatch(args),
                Cmd::KhI(args)  => khi::dispatch(args),
                Cmd::CC(args)   => cc::dispatch(args),
                Cmd::SL2(args)  => sl2::dispatch(args),
            }
        )
    }
}