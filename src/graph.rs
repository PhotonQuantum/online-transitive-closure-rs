#![allow(non_snake_case)]

use std::cell::RefCell;
use std::collections::HashSet;

use itertools::{min, Itertools};

use crate::adj_mat::AdjacencyMatrix;
use crate::error::Result;
use crate::utils::ensure_len;

#[derive(Clone, Debug, Default)]
pub struct Graph {
    M: AdjacencyMatrix,

    M_star: RefCell<AdjacencyMatrix>,

    // Condensed Graph
    V_c: HashSet<usize>,
    // V_c = {L(k)|k ∈ V}
    L: Vec<usize>,
    // L(k) = min(C(k)),                , leader of SCCs
    C: Vec<Vec<usize>>,
    // C(k) = {m ∈ V|k->*m ∧ m->*k}     , SCC, ordered C(λ)?
    P_c: Vec<HashSet<usize>>,
    // P_c(λ) = {κ ∈ V_c|(κ, λ) ∈ E_c}
    S_c: Vec<HashSet<usize>>,
    // S_c(λ) = {μ ∈ V_c|(λ, μ) ∈ E_c}
    has_cycle: bool,
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub enum EdgeEffect {
    None,
    NewEdge(Option<Vec<usize>>), // is_cycle
}

impl Graph {
    pub fn new() -> Self {
        Default::default()
    }

    pub fn with_capacity(n: usize) -> Self {
        Self {
            M: AdjacencyMatrix::with_capacity(n),
            M_star: RefCell::new(AdjacencyMatrix::with_capacity(n)),
            V_c: HashSet::with_capacity(n),
            L: Vec::with_capacity(n),
            C: Vec::with_capacity(n),
            P_c: Vec::with_capacity(n),
            S_c: Vec::with_capacity(n),
            has_cycle: false,
        }
    }

    pub fn count(&self) -> usize {
        self.M.count()
    }

    pub fn add_node(&mut self) -> usize {
        let idx = self.M.add_node();
        self.M_star.borrow_mut().add_node();
        *self.M_star.borrow_mut().get_mut(idx, idx).unwrap() = 1;
        self.V_c.insert(idx);
        ensure_len(&mut self.L, idx + 1);
        *self.L.get_mut(idx).unwrap() = idx;
        ensure_len(&mut self.C, idx + 1);
        self.C.get_mut(idx).unwrap().push(idx);
        ensure_len(&mut self.P_c, idx + 1);
        ensure_len(&mut self.S_c, idx + 1);
        idx
    }

    pub fn insert(&mut self, i: usize, j: usize) -> Result<EdgeEffect> {
        if self.M.get(i, j).unwrap() == 1 {
            return Ok(EdgeEffect::None);
        }

        *self.M.get_mut(i, j).unwrap() = 1;
        let Li = *self.L.get(i).unwrap();
        let Lj = *self.L.get(j).unwrap();

        if Li == Lj {
            Ok(EdgeEffect::NewEdge(None))
        } else {
            let new_cycle = self.M_star.borrow().get(j, i).unwrap() == 1;
            self.P_c.get_mut(Lj).unwrap().insert(Li);
            self.S_c.get_mut(Li).unwrap().insert(Lj);
            for kappa in &self.V_c {
                let kappa = *kappa;
                if self.M_star.borrow().get(kappa, i).unwrap() == 1
                    && self.M_star.borrow().get(kappa, j).unwrap() == 0
                {
                    self.adapt(kappa, Lj);
                } // R(kappa) is established
            }
            if new_cycle {
                self.has_cycle = true;
                self.join_components(j);
                Ok(EdgeEffect::NewEdge(Some(
                    self.C.get(*self.L.get(j).unwrap()).unwrap().clone(),
                )))
            } else {
                Ok(EdgeEffect::NewEdge(None))
            }
        }
    }

    fn adapt(&self, kappa: usize, first_red: usize) {
        let mut red_leaders = vec![first_red];

        while let Some(lambda) = red_leaders.pop() {
            for (k, l) in self
                .C
                .get(kappa)
                .unwrap()
                .iter()
                .cartesian_product(self.C.get(lambda).unwrap().iter())
            {
                *self.M_star.borrow_mut().get_mut(*k, *l).unwrap() = 1;
            }
            for mu in self.S_c.get(lambda).unwrap() {
                let mu = *mu;
                if self.M_star.borrow().get(kappa, mu).unwrap() == 0 {
                    red_leaders.push(mu);
                }
            }
        }
    }

    fn join_components(&mut self, j: usize) {
        let Lambda: HashSet<_> = (0..self.M.count())
            .filter(|l| {
                self.M_star.borrow().get(*l, j).unwrap() == 1
                    && self.M_star.borrow().get(j, *l).unwrap() == 1
            })
            .collect();
        let lambda = min(Lambda.clone()).unwrap();

        for kappa in &Lambda {
            let kappa = *kappa;
            *self.P_c.get_mut(kappa).unwrap() = HashSet::new();
            *self.S_c.get_mut(kappa).unwrap() = HashSet::new();
        }
        for (kappa, mu) in Lambda.iter().cartesian_product(self.V_c.iter()) {
            self.P_c.get_mut(*mu).unwrap().remove(kappa);
            self.S_c.get_mut(*mu).unwrap().remove(kappa);
        }

        for l in &Lambda {
            *self.L.get_mut(*l).unwrap() = lambda;
        }

        *self.C.get_mut(lambda).unwrap() = Lambda.iter().copied().collect();

        // adjust M_c, P_c, S_c
        for k in (0..self.M.count()).filter(|v| !Lambda.contains(v)) {
            let kappa = *self.L.get(k).unwrap();
            // k: not in Lambda, mu: in Lambda
            for mu in &Lambda {
                let mu = *mu;
                if self.M.get(k, mu).unwrap() == 1 {
                    self.S_c.get_mut(kappa).unwrap().insert(lambda);
                    self.P_c.get_mut(lambda).unwrap().insert(kappa);
                } else if self.M.get(mu, k).unwrap() == 1 {
                    self.S_c.get_mut(lambda).unwrap().insert(kappa);
                    self.P_c.get_mut(kappa).unwrap().insert(lambda);
                }
            }
        }

        // V_c := (V_c \ Lambda) \cup {lambda}
        for lam in &Lambda {
            self.V_c.remove(lam);
        }
        self.V_c.insert(lambda);
    }

    pub fn has_cycle(&self) -> bool {
        self.has_cycle
    }

    pub fn SCC(&self, reversed: bool) -> Vec<Vec<usize>> {
        let adj_list = if reversed { &self.P_c } else { &self.S_c };
        let mut V_c = self.V_c.clone();
        let mut adj_list_rev = if reversed {
            self.S_c.clone()
        } else {
            self.P_c.clone()
        };
        let mut output = vec![];
        loop {
            let in_zero_v = V_c
                .iter()
                .filter(|v| adj_list_rev.get(**v).unwrap().is_empty())
                .copied()
                .collect_vec();
            if in_zero_v.is_empty() {
                break;
            }

            for kappa in in_zero_v {
                V_c.remove(&kappa);
                output.push(self.C.get(kappa).unwrap().clone());
                for mu in adj_list.get(kappa).unwrap().iter() {
                    adj_list_rev.get_mut(*mu).unwrap().remove(&kappa);
                }
            }
        }
        output
    }
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use crate::graph::{EdgeEffect, Graph};

    #[test]
    fn adj_graph_works() {
        println!("adj graph works");
        let mut g = Graph::new();
        for i in 0..10 {
            assert_eq!(g.add_node(), i);
            assert_eq!(g.M_star.borrow().get(i, i).unwrap(), 1);
        }
        assert_eq!(g.insert(1, 2).unwrap(), EdgeEffect::NewEdge(None));
        assert_eq!(g.insert(2, 4).unwrap(), EdgeEffect::NewEdge(None));
        assert_eq!(g.insert(4, 3).unwrap(), EdgeEffect::NewEdge(None));
        let r = g.insert(3, 1).unwrap();
        if let EdgeEffect::NewEdge(Some(l)) = r {
            assert_eq!(l.iter().sorted().collect_vec(), [&1, &2, &3, &4]);
        } else {
            assert!(false);
        }
        assert_eq!(g.insert(5, 6).unwrap(), EdgeEffect::NewEdge(None));
        assert_eq!(g.insert(6, 7).unwrap(), EdgeEffect::NewEdge(None));
        let r = g.insert(7, 5).unwrap();
        if let EdgeEffect::NewEdge(Some(l)) = r {
            assert_eq!(l.iter().sorted().collect_vec(), [&5, &6, &7]);
        } else {
            assert!(false);
        }
        assert_eq!(g.insert(4, 5).unwrap(), EdgeEffect::NewEdge(None));
        assert_eq!(g.insert(2, 3).unwrap(), EdgeEffect::NewEdge(None));
        assert_eq!(g.insert(7, 8).unwrap(), EdgeEffect::NewEdge(None));
        assert_eq!(g.insert(8, 9).unwrap(), EdgeEffect::NewEdge(None));
        let r = g.insert(9, 7).unwrap();
        if let EdgeEffect::NewEdge(Some(l)) = r {
            assert_eq!(l.iter().sorted().collect_vec(), [&5, &6, &7, &8, &9]);
        } else {
            assert!(false);
        }
        assert_eq!(g.insert(7, 0).unwrap(), EdgeEffect::NewEdge(None));
        assert_eq!(g.insert(4, 0).unwrap(), EdgeEffect::NewEdge(None));
        assert_eq!(g.insert(4, 3).unwrap(), EdgeEffect::None);

        assert_eq!(g.V_c.iter().sorted().collect_vec(), [&0, &1, &5]);
        assert_eq!(
            g.V_c
                .iter()
                .sorted()
                .map(|v| (v, g.C.get(*v).unwrap().iter().sorted().collect_vec()))
                .collect_vec(),
            vec![
                (&0, vec![&0]),
                (&1, vec![&1, &2, &3, &4]),
                (&5, vec![&5, &6, &7, &8, &9])
            ]
        );
        assert_eq!(
            g.V_c
                .iter()
                .sorted()
                .map(|v| (v, g.S_c.get(*v).unwrap().iter().sorted().collect_vec()))
                .collect_vec(),
            vec![(&0, vec![]), (&1, vec![&0, &5]), (&5, vec![&0])]
        );
        assert_eq!(
            g.SCC(false)
                .iter()
                .map(|s| s.iter().sorted().collect_vec())
                .collect_vec(),
            vec![vec![&1, &2, &3, &4], vec![&5, &6, &7, &8, &9], vec![&0]]
        );
        assert_eq!(
            g.SCC(true)
                .iter()
                .map(|s| s.iter().sorted().collect_vec())
                .collect_vec(),
            vec![vec![&0], vec![&5, &6, &7, &8, &9], vec![&1, &2, &3, &4]]
        );
        assert_eq!(
            g.M.edges().collect_vec(),
            vec![
                (1, 2),
                (2, 3),
                (2, 4),
                (3, 1),
                (4, 0),
                (4, 3),
                (4, 5),
                (5, 6),
                (6, 7),
                (7, 0),
                (7, 5),
                (7, 8),
                (8, 9),
                (9, 7)
            ]
        );
    }
}
