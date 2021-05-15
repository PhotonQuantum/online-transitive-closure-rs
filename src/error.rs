use thiserror::Error;

pub(crate) type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, Error, Copy, Clone, Eq, PartialEq)]
pub enum Error {
    #[error("Invalid index")]
    InvalidIndex,
}
