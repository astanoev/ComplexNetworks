package cs.graph;

import java.awt.Color;
import java.util.ArrayList;
import java.util.HashMap;

public class MultiNode extends Node {
	HashMap<String, Node> nodes = new HashMap<String, Node>();
	public Color color;
	
	public float[] getColor(){
		return new float[]{color.getRed()/256.0f, color.getGreen()/256.0f,color.getBlue()/256.0f};
	}
		
	public MultiNode(int nodeID) {
		super(nodeID);
		color = Color.getHSBColor(0.6f + (0.618033988749895f * id) % 1, .75f, 0.75f);
	}
	
	public void addNode(Node n, String type){
		nodes.put(type, n);
		n.putAttr("superNode", this);
	}
	
	public Node getNode(String type){
		if (nodes.get(type)!=null)
			return nodes.get(type);
		else
			return new Node(-1);
	}
	
	public ArrayList<Edge> outEdges(String type){
		Node n = nodes.get(type);
		if(n==null)
			return new ArrayList<Edge>();
		return n.outEdges;
	}
}
