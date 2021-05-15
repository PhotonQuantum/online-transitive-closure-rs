mod adj_mat;
mod error;
mod graph;
mod utils;

pub use adj_mat::AdjacencyMatrix;
pub use graph::{Graph, EdgeEffect};
pub use error::Error;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
