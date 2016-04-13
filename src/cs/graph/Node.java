package cs.graph;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;

public class Node implements Comparable<Integer> {
	public ArrayList<Edge> inEdges = new ArrayList<Edge>();
	public ArrayList<Edge> outEdges = new ArrayList<Edge>();
	private HashMap<String, Object> attr = null;
	public int id = -1;

	public Node(int id) {
		this.id = id;
	}
	
	public Object getAttr(String key) {
		if (attr == null) {
			attr = new HashMap<String, Object>();
		}
		return attr.get(key);
	}

	public boolean hasAttr(String key) {
		if (attr == null) {
			attr = new HashMap<String, Object>();
		}
		return attr.containsKey(key);
	}

	public void putAttr(String key, Object value) {
		if (attr == null) {
			attr = new HashMap<String, Object>();
		}
		attr.put(key, value);
	}

	protected Edge addEdge(Node dest) {
		Edge edge = new Edge(this, dest);
		int index = Collections.binarySearch(outEdges, edge);
		if (index >= 0) {
			// System.out.println("The Edge (" + this.id + "," + dest.id
			// + ") already exists");
			return null;//outEdges.get(index);
		} else {
			outEdges.add(-index - 1, edge);
			return edge;
		}
	}

	public boolean hasEdge(Node dest) {
		Edge edge = new Edge(this, dest);
		int index = Collections.binarySearch(outEdges, edge);
		if (index >= 0) {
			return true;
		} else {
			return false;
		}
	}

	public Edge getEdge(Node dest) {
		Edge edge = new Edge(this, dest);
		int index = Collections.binarySearch(outEdges, edge);
		if (index >= 0) {
			return outEdges.get(index);
		} else {
			return null;
		}
	}

	private Comparator<Edge> comparator = new Comparator<Edge>() {
		@Override
		public int compare(Edge o1, Edge o2) {
			return o1.src.id - o2.src.id;
		}
	};

	protected void addInEdge(Edge edge) {
		int index = Collections.binarySearch(inEdges, edge, comparator);
		inEdges.add(-index - 1, edge);
	}

	@Override
	public int compareTo(Integer nodeid) {
		return this.id - nodeid;
	}

	@Override
	public int hashCode() {
		return this.id;
	}

	@Override
	public boolean equals(Object obj) {
		if (obj instanceof Integer) {
			return ((Integer) obj) == this.id;
		} else if (obj instanceof Node) {
			return ((Node) obj).id == this.id;
		} else {
			return super.equals(obj);
		}
	}

	@Override
	public String toString() {
		if (this.hasAttr("label")) {
			return (String) this.getAttr("label");
		} else
			return this.id + "";
	}
}
