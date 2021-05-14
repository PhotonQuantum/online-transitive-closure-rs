use thiserror::Error;

pub(crate) type Result<T> = std::result::Result<T, GraphError>;

#[derive(Debug, Error)]
pub enum GraphError {
    #[error("Invalid index")]
    InvalidIndex,
}
