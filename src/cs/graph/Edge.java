package cs.graph;

import java.util.HashMap;

public class Edge implements Comparable<Edge> {
	private HashMap<String, Object> attr = null;
	public Node dest;
	public Node src;

	public Object getAttr(String key) {
		if (attr == null) {
			attr = new HashMap<String, Object>();
		}
		return attr.get(key);
	}

	public Object getAttr(String key, Object startValue) {
		if (attr == null) {
			attr = new HashMap<String, Object>();
		}
		if (attr.get(key) == null) {
			attr.put(key, startValue);
		}
		return attr.get(key);
	}

	public void putAttr(String key, Object value) {
		if (attr == null) {
			attr = new HashMap<String, Object>();
		}
		attr.put(key, value);
	}

	public Edge(Node src, Node dest) {
		this.src = src;
		this.dest = dest;
	}

	@Override
	public int compareTo(Edge edge) {
		return this.dest.id - edge.dest.id;
	}
}
