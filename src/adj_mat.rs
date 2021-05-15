use std::cmp;
use std::iter::Enumerate;
use std::slice::Iter;

use crate::error::*;
use crate::utils::ensure_len;

#[derive(Debug, Clone)]
pub struct AdjacencyMatrix {
    mat: Vec<usize>,
    nodes_count: usize,
    capacity: usize,
}

#[derive(Debug, Clone)]
pub struct Edges<'a> {
    source: Enumerate<Iter<'a, usize>>,
    capacity: usize
}

impl<'a> From<&'a AdjacencyMatrix> for Edges<'a> {
    fn from(mat: &'a AdjacencyMatrix) -> Self {
        Self {
            source: mat.mat.iter().enumerate(),
            capacity: mat.capacity
        }
    }
}

impl Iterator for Edges<'_> {
    type Item = (usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        self.source
            .find(|(_, weight)| **weight > 0)
            .map(|(idx, _)| AdjacencyMatrix::to_pos(idx, self.capacity))
    }
}

impl Default for AdjacencyMatrix {
    fn default() -> Self {
        let mut mat = Self {
            mat: vec![],
            nodes_count: 0,
            capacity: 0,
        };
        mat.extend_matrix(0, false);
        mat
    }
}

impl AdjacencyMatrix {
    pub fn new() -> Self {
        Default::default()
    }

    pub fn count(&self) -> usize {
        self.nodes_count
    }

    pub fn with_capacity(capacity: usize) -> Self {
        let mut mat = Self {
            mat: vec![],
            nodes_count: 0,
            capacity: 0,
        };
        mat.extend_matrix(capacity, true);
        mat
    }

    pub fn add_node(&mut self) -> usize {
        self.nodes_count += 1;
        if self.nodes_count > self.capacity {
            self.extend_matrix(self.nodes_count, false);
        }
        self.nodes_count - 1
    }

    pub fn ensure_nodes(&mut self, n: usize) {
        if self.nodes_count < n {
            self.nodes_count = n;
            self.extend_matrix(self.nodes_count, false);
        }
    }

    pub fn get(&self, i: usize, j: usize) -> Result<usize> {
        if i >= self.nodes_count || j >= self.nodes_count {
            return Err(GraphError::InvalidIndex);
        }
        let pos = self.to_linear_pos(i, j);
        Ok(*self.mat.get(pos).unwrap())
    }

    pub unsafe fn get_unchecked(&self, i: usize, j: usize) -> usize {
        let pos = self.to_linear_pos(i, j);
        *self.mat.get_unchecked(pos)
    }

    pub fn get_mut(&mut self, i: usize, j: usize) -> Result<&mut usize> {
        if i >= self.nodes_count || j >= self.nodes_count {
            return Err(GraphError::InvalidIndex);
        }
        let pos = self.to_linear_pos(i, j);
        Ok(self.mat.get_mut(pos).unwrap())
    }

    pub unsafe fn get_unchecked_mut(&mut self, i: usize, j: usize) -> &mut usize {
        let pos = self.to_linear_pos(i, j);
        self.mat.get_unchecked_mut(pos)
    }

    pub fn edges(&self) -> Edges {
        Edges::from(self)
    }

    fn extend_matrix(&mut self, new_capacity: usize, exact: bool) {
        self.capacity =
            Self::extend_flat_square_matrix(&mut self.mat, self.capacity, new_capacity, exact);
    }

    // adapted from petgraph
    #[inline]
    fn extend_flat_square_matrix<T: Default>(
        node_adjacencies: &mut Vec<T>,
        old_node_capacity: usize,
        new_node_capacity: usize,
        exact: bool,
    ) -> usize {
        // Grow the capacity by exponential steps to avoid repeated allocations.
        // Disabled for the with_capacity constructor.
        let new_node_capacity = if exact {
            new_node_capacity
        } else {
            const MIN_CAPACITY: usize = 4;
            cmp::max(new_node_capacity.next_power_of_two(), MIN_CAPACITY)
        };

        // Optimization: when resizing the matrix this way we skip the first few grows to make
        // small matrices a bit faster to work with.

        ensure_len(node_adjacencies, new_node_capacity.pow(2));
        for c in (1..old_node_capacity).rev() {
            let pos = c * old_node_capacity;
            let new_pos = c * new_node_capacity;

            // Move the slices directly if they do not overlap with their new position
            if pos + old_node_capacity <= new_pos {
                let old = node_adjacencies[pos..pos + old_node_capacity].as_mut_ptr();
                let new = node_adjacencies[new_pos..new_pos + old_node_capacity].as_mut_ptr();

                // SAFE: new starts at least `old_node_capacity` positions after old.
                unsafe {
                    std::ptr::swap_nonoverlapping(old, new, old_node_capacity);
                }
            } else {
                for i in (0..old_node_capacity).rev() {
                    node_adjacencies.as_mut_slice().swap(pos + i, new_pos + i);
                }
            }
        }

        new_node_capacity
    }

    fn to_linear_pos(&self, row: usize, col: usize) -> usize {
        row * self.capacity + col
    }

    fn to_pos(idx: usize, width: usize) -> (usize, usize) {
        (idx / width, idx % width)
    }
}

#[cfg(test)]
mod tests {
    use crate::adj_mat::AdjacencyMatrix;

    #[test]
    fn adj_mat_works() {
        let mut mat = AdjacencyMatrix::new();
        for i in 0..10 {
            assert_eq!(mat.add_node(), i);
        }
        assert_eq!(mat.get(0, 0).unwrap(), 0);
        assert_eq!(mat.get(3, 1).unwrap(), 0);
        assert_eq!(mat.get(3, 9).unwrap(), 0);
        *mat.get_mut(0, 0).unwrap() = 1;
        *mat.get_mut(3, 1).unwrap() = 1;
        assert_eq!(mat.get(0, 0).unwrap(), 1);
        assert_eq!(mat.get(3, 1).unwrap(), 1);
        *mat.get_mut(5, 1).unwrap() = 1;
        *mat.get_mut(9, 7).unwrap() = 1;
        assert_eq!(mat.get(5, 1).unwrap(), 1);
        assert_eq!(mat.get(9, 7).unwrap(), 1);
        *mat.get_mut(5, 1).unwrap() = 0;
        assert_eq!(mat.get(5, 1).unwrap(), 0);
        assert!(mat.get_mut(0, 10).is_err());
        assert!(mat.get(0, 10).is_err());
    }

    #[test]
    fn adj_mat_capacity() {
        let mat = AdjacencyMatrix::with_capacity(1);
        assert_eq!(mat.capacity, 1);
        let mut mat = AdjacencyMatrix::with_capacity(5);
        assert_eq!(mat.capacity, 5);
        for _ in 0..6 {
            mat.add_node();
        }
        assert_eq!(mat.capacity, 8);
    }
}
