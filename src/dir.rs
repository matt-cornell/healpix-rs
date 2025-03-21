//! Direction enumerations.

/// A positive or negative value
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(i8)]
pub enum PosOrNeg {
    Neg = -1,
    Zero = 0,
    Pos = 1,
}

/// A cardinal direction.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum Cardinal {
    N,
    E,
    S,
    W,
}
impl Cardinal {
    pub const VALUES: [Self; 4] = [Self::N, Self::E, Self::S, Self::W];
    /// Convert this to an index. This starts with 0 as north and rotates clockwise.
    #[inline(always)]
    pub const fn index(self) -> u8 {
        self as u8
    }
    /// Convert an index to a cardinal direction, using the ordering specified in [`index`].
    #[inline(always)]
    pub const fn from_index(index: u8) -> Self {
        match index {
            0 => Self::N,
            1 => Self::E,
            2 => Self::S,
            3 => Self::W,
            _ => panic!("Cardinal index out of bounds"),
        }
    }
    /// Attempt to convert an index to a cardinal direction, using the ordering specified in [`index`].
    #[inline(always)]
    pub const fn try_from_index(index: u8) -> Option<Self> {
        match index {
            0 => Some(Self::N),
            1 => Some(Self::E),
            2 => Some(Self::S),
            3 => Some(Self::W),
            _ => None,
        }
    }
    /// Get the direction of the X direction. East is positive.
    #[inline(always)]
    pub const fn x(self) -> PosOrNeg {
        match self {
            Self::W => PosOrNeg::Neg,
            Self::E => PosOrNeg::Pos,
            _ => PosOrNeg::Zero,
        }
    }
    /// Get the direction of the Y direction. North is positive.
    #[inline(always)]
    pub const fn y(self) -> PosOrNeg {
        match self {
            Self::S => PosOrNeg::Neg,
            Self::N => PosOrNeg::Pos,
            _ => PosOrNeg::Zero,
        }
    }
    /// Get the x and y directions of this coordinate. East is positive X and north is positive Y.
    #[inline(always)]
    pub const fn xy(self) -> [PosOrNeg; 2] {
        match self {
            Self::N => [PosOrNeg::Zero, PosOrNeg::Pos],
            Self::E => [PosOrNeg::Pos, PosOrNeg::Zero],
            Self::S => [PosOrNeg::Zero, PosOrNeg::Neg],
            Self::W => [PosOrNeg::Neg, PosOrNeg::Zero],
        }
    }
    /// Try to create find a cardinal direction from X- and Y- directions.
    #[inline(always)]
    pub const fn from_xy(x: PosOrNeg, y: PosOrNeg) -> Option<Self> {
        match [x, y] {
            [PosOrNeg::Zero, PosOrNeg::Pos] => Some(Self::N),
            [PosOrNeg::Pos, PosOrNeg::Zero] => Some(Self::E),
            [PosOrNeg::Zero, PosOrNeg::Neg] => Some(Self::S),
            [PosOrNeg::Neg, PosOrNeg::Zero] => Some(Self::W),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum Ordinal {
    NE,
    SE,
    SW,
    NW,
}
impl Ordinal {
    pub const VALUES: [Self; 4] = [Self::NE, Self::SE, Self::SW, Self::NW];
    /// Convert this to an index. This starts with 0 as northeast and rotates clockwise.
    /// The value is 45 degrees to the right of the cardinal value with the same index.
    #[inline(always)]
    pub const fn index(self) -> u8 {
        self as u8
    }
    /// Convert an index to a cardinal direction, using the ordering specified in [`index`].
    #[inline(always)]
    pub const fn from_index(index: u8) -> Self {
        match index {
            0 => Self::NE,
            1 => Self::SE,
            2 => Self::SW,
            3 => Self::NW,
            _ => panic!("Ordinal index out of bounds"),
        }
    }
    /// Attempt to convert an index to a cardinal direction, using the ordering specified in [`index`].
    #[inline(always)]
    pub const fn try_from_index(index: u8) -> Option<Self> {
        match index {
            0 => Some(Self::NE),
            1 => Some(Self::SE),
            2 => Some(Self::SW),
            3 => Some(Self::NW),
            _ => None,
        }
    }
    /// Get the direction of the X direction. East is positive.
    #[inline(always)]
    pub const fn x(self) -> PosOrNeg {
        match self {
            Self::NW | Self::SW => PosOrNeg::Neg,
            Self::NE | Self::SE => PosOrNeg::Pos,
        }
    }
    /// Get the direction of the Y direction. North is positive.
    #[inline(always)]
    pub const fn y(self) -> PosOrNeg {
        match self {
            Self::SW | Self::SE => PosOrNeg::Neg,
            Self::NW | Self::NE => PosOrNeg::Pos,
        }
    }
    /// Get the x and y directions of this coordinate. East is positive X and north is positive Y.
    #[inline(always)]
    pub const fn xy(self) -> [PosOrNeg; 2] {
        match self {
            Self::NE => [PosOrNeg::Pos, PosOrNeg::Pos],
            Self::SE => [PosOrNeg::Pos, PosOrNeg::Neg],
            Self::SW => [PosOrNeg::Neg, PosOrNeg::Neg],
            Self::NW => [PosOrNeg::Neg, PosOrNeg::Pos],
        }
    }
    /// Try to create find a ordinal direction from X- and Y- directions.
    #[inline(always)]
    pub const fn from_xy(x: PosOrNeg, y: PosOrNeg) -> Option<Self> {
        match [x, y] {
            [PosOrNeg::Pos, PosOrNeg::Pos] => Some(Self::NE),
            [PosOrNeg::Pos, PosOrNeg::Neg] => Some(Self::SE),
            [PosOrNeg::Neg, PosOrNeg::Neg] => Some(Self::SW),
            [PosOrNeg::Neg, PosOrNeg::Pos] => Some(Self::NW),
            _ => None,
        }
    }

    /// Rotate clockwise by a given amount.
    #[inline(always)]
    pub const fn rot_cw(self, dist: u8) -> Self {
        Self::from_index(self.index().wrapping_add(dist) % 4)
    }
    /// Rotate counter-clockwise by a given amount.
    #[inline(always)]
    pub const fn rot_ccw(self, dist: u8) -> Self {
        Self::from_index(self.index().wrapping_sub(dist) % 4)
    }
    /// Clockwise rotation distance to rotate self to the other direction.
    #[inline(always)]
    pub const fn dist_cw(self, other: Self) -> u8 {
        other.index().wrapping_sub(self.index()) % 4
    }
    /// Counter-clockwise rotation distance to rotate self to the other direction.
    #[inline(always)]
    pub const fn dist_ccw(self, other: Self) -> u8 {
        self.index().wrapping_sub(other.index()) % 4
    }
}

/// A direction, directed
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum Direction {
    N,
    NE,
    E,
    SE,
    S,
    SW,
    W,
    NW,
}
impl Direction {
    pub const VALUES: [Self; 8] = [
        Self::N,
        Self::NE,
        Self::E,
        Self::SE,
        Self::S,
        Self::SW,
        Self::W,
        Self::NW,
    ];
    #[inline(always)]
    pub const fn index(self) -> u8 {
        self as u8
    }
    #[inline(always)]
    pub const fn from_index(index: u8) -> Self {
        match index {
            0 => Self::N,
            1 => Self::NE,
            2 => Self::E,
            3 => Self::SE,
            4 => Self::S,
            5 => Self::SW,
            6 => Self::W,
            7 => Self::NW,
            _ => panic!("Direction index out of bounds"),
        }
    }
    #[inline(always)]
    pub const fn try_from_index(index: u8) -> Option<Self> {
        match index {
            0 => Some(Self::N),
            1 => Some(Self::NE),
            2 => Some(Self::E),
            3 => Some(Self::SE),
            4 => Some(Self::S),
            5 => Some(Self::SW),
            6 => Some(Self::W),
            7 => Some(Self::NW),
            _ => None,
        }
    }
    /// Get the direction of the X direction. East is positive.
    #[inline(always)]
    pub const fn x(self) -> PosOrNeg {
        match self {
            Self::NW | Self::SW | Self::W => PosOrNeg::Neg,
            Self::NE | Self::SE | Self::E => PosOrNeg::Pos,
            Self::N | Self::S => PosOrNeg::Zero,
        }
    }
    /// Get the direction of the Y direction. North is positive.
    #[inline(always)]
    pub const fn y(self) -> PosOrNeg {
        match self {
            Self::SW | Self::SE | Self::S => PosOrNeg::Neg,
            Self::NW | Self::NE | Self::N => PosOrNeg::Pos,
            Self::W | Self::E => PosOrNeg::Zero,
        }
    }
    /// Get the x and y directions of this coordinate. East is positive X and north is positive Y.
    #[inline(always)]
    pub const fn xy(self) -> [PosOrNeg; 2] {
        match self {
            Self::N => [PosOrNeg::Zero, PosOrNeg::Pos],
            Self::NE => [PosOrNeg::Pos, PosOrNeg::Pos],
            Self::E => [PosOrNeg::Pos, PosOrNeg::Zero],
            Self::SE => [PosOrNeg::Pos, PosOrNeg::Neg],
            Self::S => [PosOrNeg::Zero, PosOrNeg::Neg],
            Self::SW => [PosOrNeg::Neg, PosOrNeg::Neg],
            Self::W => [PosOrNeg::Neg, PosOrNeg::Zero],
            Self::NW => [PosOrNeg::Neg, PosOrNeg::Pos],
        }
    }
    /// Try to create find a cardinal direction from X- and Y- directions.
    #[inline(always)]
    pub const fn from_xy(x: PosOrNeg, y: PosOrNeg) -> Option<Self> {
        match [x, y] {
            [PosOrNeg::Zero, PosOrNeg::Pos] => Some(Self::N),
            [PosOrNeg::Pos, PosOrNeg::Pos] => Some(Self::NE),
            [PosOrNeg::Pos, PosOrNeg::Zero] => Some(Self::E),
            [PosOrNeg::Pos, PosOrNeg::Neg] => Some(Self::SE),
            [PosOrNeg::Zero, PosOrNeg::Neg] => Some(Self::S),
            [PosOrNeg::Neg, PosOrNeg::Neg] => Some(Self::SW),
            [PosOrNeg::Neg, PosOrNeg::Zero] => Some(Self::W),
            [PosOrNeg::Neg, PosOrNeg::Pos] => Some(Self::NW),
            _ => None,
        }
    }
}

impl From<Cardinal> for Direction {
    fn from(value: Cardinal) -> Self {
        match value {
            Cardinal::N => Self::N,
            Cardinal::E => Self::E,
            Cardinal::S => Self::S,
            Cardinal::W => Self::W,
        }
    }
}

impl From<Ordinal> for Direction {
    fn from(value: Ordinal) -> Self {
        match value {
            Ordinal::NE => Self::NE,
            Ordinal::SE => Self::SE,
            Ordinal::SW => Self::SW,
            Ordinal::NW => Self::NW,
        }
    }
}
