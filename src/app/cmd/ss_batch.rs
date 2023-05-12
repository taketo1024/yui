#![allow(unused)]

use std::fmt::Display;

use indexmap::IndexMap;
use itertools::Itertools;
use log::error;
use serde::Deserialize;
use csv::Writer;
use yui_link::{Edge, Link};
use super::ss::{Args as SSArgs, run as ss_run};
use crate::utils::*;

#[derive(Debug, clap::Args)]
pub struct Args { 
    input: String,
    output: String,

    #[arg(long)]
    pub debug: bool
}

type PDCode = Vec<[Edge; 4]>;

#[derive(Deserialize, Debug)]
struct InputData { 
    configs: Vec<Config>,
    targets: IndexMap<String, String>
}

#[derive(Deserialize, Debug)]
struct Config {
    name: String,
    c_value: String,
    c_type: CType,

    #[serde(default)]
    reduced: bool,

    #[serde(default)]
    mirror: bool,

    #[serde(default)]
    order: Option<usize>,
}

pub fn run(args: &Args) -> Result<String, Box<dyn std::error::Error>> {
    let json_str = std::fs::read_to_string(&args.input)?;
    let InputData{ configs, targets } = serde_json::from_str(&json_str)?;
    
    let mut wtr = csv_writer(&args.output)?;

    let head = [
        vec!["name".to_string()], 
        configs.iter().map(|c| c.name.to_string()).collect()
    ].concat();

    wtr.write_record(head);
    wtr.flush()?;
    
    for (name, link) in targets.iter() { 
        let mut row = vec![name.to_string()];

        for config in configs.iter() { 
            let ss_args = make_args(&link, &config, args.debug);
            let ss_res = ss_run(&ss_args);

            let val = if let Ok(val) = ss_res { 
                val.to_string()
            } else { 
                if let Err(e) = ss_res {
                    error!("{}", e);
                    eprintln!("{e}");
                }
                "!".to_string()            
            };

            row.push(val);
        }
        
        wtr.write_record(row);
        wtr.flush()?;    
    }
    
    let msg = format!("result saved to: {}", &args.output);
    Ok(msg)
}

fn make_args(link: &String, config: &Config, debug: bool) -> SSArgs { 
    SSArgs {
        link:    link.clone(),
        c_value: config.c_value.clone(),
        c_type:  config.c_type.clone(),
        mirror:  config.mirror,
        order:   config.order.unwrap_or(1),
        reduced: config.reduced,
        old:     false,
        debug:   debug,
    }
}