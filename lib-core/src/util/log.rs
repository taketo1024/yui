pub fn init_simple_logger(l: log::LevelFilter) -> Result<(), log::SetLoggerError> { 
    use simplelog::*;

    let mut cb = simplelog::ConfigBuilder::new();
    cb.set_location_level(LevelFilter::Off);
    cb.set_target_level(LevelFilter::Off);
    cb.set_thread_level(LevelFilter::Off);
    cb.set_level_color(Level::Trace, Some(Color::Green));
    let config = cb.build();

    TermLogger::init(
        l,
        config,
        TerminalMode::Mixed,
        ColorChoice::Always
    )
}

