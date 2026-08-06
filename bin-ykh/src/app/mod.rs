//! The CLI application: argument parsing, command dispatch and the commands.

mod app;
pub use app::{App, CliArgs};

mod args;
mod err;
mod cmd;
mod utils;