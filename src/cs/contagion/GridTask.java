package cs.contagion;

import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.Callable;

import cs.graph.Edge;
import cs.graph.Graph;
import cs.graph.MultiGraph;
import cs.graph.MultiNode;
import cs.graph.Node;

public class GridTask implements Callable<List<Node>> {
	private Graph g;
	private int i_start;
	private int i_end;
	private MultiGraph mg;

	public GridTask() {
	}

	public GridTask(Graph g, int i_start, int i_end, MultiGraph mg) {
		this.g = g;
		this.i_start = i_start;
		this.i_end = i_end;
		this.mg = mg;
	}

	public List<Node> call() {

		Random ran = new SecureRandom();
		for (Node n : g.nodes.subList(i_start, i_end)) {
			try {
				MultiNode mn = (MultiNode) n.getAttr("curCol");
				n.putAttr("nextCol", mn);
				double probTotal = 0;
				double probGen = ran.nextDouble();
				for (Edge e : mn.outEdges("delta")) {
					double edgeValue = ((Double) e.getAttr("value"));
					if (probGen <= probTotal + edgeValue) {
						n.putAttr("nextCol", mg.getAddNode(e.dest.id));
						probTotal += edgeValue;
						break;
					}
					probTotal += edgeValue;
				}
				if (probGen <= probTotal)// flag
					continue;

				ArrayList<Node> infectiveNeighbors = new ArrayList<Node>();
				for (Edge e : n.outEdges) {
					Edge ed = mn.getNode("beta").getEdge(((MultiNode) e.dest.getAttr("curCol")).getNode("beta"));
					if (ed != null)
						if (ran.nextDouble() <= (Double) ed.getAttr("value"))
							infectiveNeighbors.add(e.dest);
				}
				if (infectiveNeighbors.size() > 0)
					n.putAttr("nextCol", (MultiNode) infectiveNeighbors.get(ran.nextInt(infectiveNeighbors.size())).getAttr("curCol"));
			} catch (Exception ex) {
				System.out.println("error");
			}
		}

		return g.nodes.subList(i_start, i_end);
	}

}

