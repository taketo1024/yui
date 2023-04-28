#[derive(PartialEq, Eq, Clone)]
pub enum Color { A, B }

impl Color { 
    pub fn is_a(&self) -> bool { 
        self == &Color::A
    }

    pub fn other(&self) -> Self { 
        match self { 
            Color::A => Color::B,
            Color::B => Color::A
        }
    }
}