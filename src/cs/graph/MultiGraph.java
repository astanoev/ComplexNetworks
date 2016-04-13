package cs.graph;

import java.awt.Color;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Random;

public class MultiGraph {

	public MultiGraph(){};
	
	public MultiGraph(MultiGraph other) {
		this.nedges = other.nedges;
		this.nodes = other.nodes;
		this.attr = other.attr;
	}

	public ArrayList<MultiNode> nodes = new ArrayList<MultiNode>();
	public ArrayList<Graph> graphs = new ArrayList<Graph>();

	public int nedges = 0;

	public static final boolean PRINT_PROGRESS = true;

	private HashMap<String, Object> attr = new HashMap<String, Object>();

	public void addGraph(Graph g){
		String type = (String)g.getAttr("type");
		for(Node n: g.nodes){
			MultiNode superNode = getAddNode(n.id);
			superNode.addNode(n, type);
		}
	}
	
	public Graph getGraph(String type){
		for(Graph g: graphs){
			if(((String)g.getAttr("type")).equalsIgnoreCase(type))
				return g;
		}
		return null;
	}
	
	public MultiNode getAddNode(int nodeID) {
		int index = Collections.binarySearch(nodes, nodeID);
		if (index >= 0) {
			return nodes.get(index);
		} else {
			MultiNode tmp = new MultiNode(nodeID);
			nodes.add(-index - 1, tmp);
			return tmp;
		}
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
	
	public static MultiGraph createFullMeshDiseaseGraph(int M, int N,
			double beta_mean, double sd) {
		Random ran = new SecureRandom();
		MultiGraph mg = new MultiGraph();
		Graph g = new Graph(false);
		g.putAttr("type", "beta");
		for (int i = 0; i < M; i++) {
			Node n = new Node(i);
			g.nodes.add(n);
			for (Node n2 : g.nodes) {
				Edge e = g.addEdge(n2.id, n.id);
				e.putAttr("value", beta_mean + sd * ran.nextGaussian());
				if (n.id != n2.id) {
					Edge e2 = g.addEdge(n.id, n2.id);
					e2.putAttr("value", beta_mean + sd * ran.nextGaussian());
				}
			}
		}
		mg.addGraph(g);
		return mg;
	}
	
	public static MultiGraph FullMeshDiseaseToyGraph(int M, int N) {
		//double[][] beta = {{0.1797, 0.1807, 0.0964, 0.15}, {0.1357, 0.0963, 0.2939, 0.02}, {0.2892, 0.1271, 0.0989, 0.17}, {0.1, 0.1, 0.04, 0.12}};
		//double[][] beta = {{0.01797, 0.01807, 0.00964, 0.015}, {0.01357, 0.00963, 0.02939, 0.002}, {0.02892, 0.01271, 0.00989, 0.017}, {0.01, 0.01, 0.004, 0.012}};
		//double[][] beta = {{0.1797, 0.1807, 0.0964}, {0.1357, 0.0963, 0.2939}, {0.2892, 0.1271, 0.0989}};
		//double[][] beta = {{0.18, 0.18, 0.095}, {0.135, 0.095, 0.295}, {0.29, 0.125, 0.1}};
		double[][] beta = {{0.225, 0.225, 0.123}, {0.17, 0.145, 0.363}, {0.362, 0.156, 0.125}};
		//double[][] beta = {{0.025, 0.09, 0.025}, {0.025, 0.025, 0.05}, {0.15, 0.025, 0.025}};***
		//double[][] beta = {{0.125, 0.45, 0.125}, {0.125, 0.125, 0.25}, {0.75, 0.125, 0.125}};
		//double[][] beta = {{0, 0.05, 0.07}, {0.07, 0, 0.05}, {0.05, 0.07, 0}};
		//double[][] beta = {{0, 0.015, 0.17}, {0.07, 0.11, 0.05}, {0.035, 0.1, 0.002}};
		MultiGraph mg = new MultiGraph();
		Graph g = new Graph(false);
		g.putAttr("type", "beta");
		for (int i = 0; i < M; i++) {
			Node n = new Node(i);
			g.nodes.add(n);
			for (Node n2 : g.nodes) {
				if(beta[n2.id][n.id]>0){
					Edge e = g.addEdge(n2.id, n.id);
					e.putAttr("value", beta[n2.id][n.id]);
				}
				if (n.id != n2.id) {
					if(beta[n.id][n2.id]>0){
						Edge e2 = g.addEdge(n.id, n2.id);					
						e2.putAttr("value", beta[n.id][n2.id]);
					}
				}
			}
		}
		mg.addGraph(g);
		return mg;
	}

	public static MultiGraph LanguageSpread(int N) {
		//double[][] beta = {{0.1797, 0.1807, 0.0964, 0.15}, {0.1357, 0.0963, 0.2939, 0.02}, {0.2892, 0.1271, 0.0989, 0.17}, {0.1, 0.1, 0.04, 0.12}};
		double[][] beta = {{0.18, 0.18, 0.1, 0.15}, {0.14, 0.1, 0.3, 0.02}, {0.29, 0.13, 0.1, 0.17}, {0.1, 0.1, 0.04, 0.12}};
		double betaO = 1;
		int M = beta.length;
		int L = (int)Math.sqrt(N);
		//double[][] beta = {{0.25, 0.25}, {0.25, 0.25}};
		//double[][] beta = {{0.9, 0.1}, {0.1, 0.9}};
		int[][] sources = new int[][] { { L / 4, L / 4 },
				{ 3 * L / 4, 3 * L / 4 }, { L / 4, 3 * L / 4 },
				{ 3 * L / 4, L / 4 } };
		MultiGraph mg = new MultiGraph();
		Graph g = new Graph(false);
		for (int i = 0; i < M; i++) {
			Node n = new Node(i);
			g.nodes.add(n);
			n.putAttr("sourceCoords", sources[i]);
			for (Node n2 : g.nodes) {
				Edge e = g.addEdge(n2.id, n.id);
				e.putAttr("Type", 0);
				e.putAttr("Value", betaO*beta[n2.id][n.id]);
				if (n.id != n2.id) {
					Edge e2 = g.addEdge(n.id, n2.id);
					e2.putAttr("Type", 0);
					e2.putAttr("Value", betaO*beta[n.id][n2.id]);
				}
				else{
					e.putAttr("Type", 1);
					e.putAttr("Value", 1.0);
				}
			}
		}
		//double betaSn = 0.5;
//		Node sn = new Node(M);
//		Col colSn = new Col(M, 0, N, 0);
//		colSn.setColor(new float[] { 1, 1, 1 });
//		g.nodes.add(sn);
//		sn.putAttr("Col",colSn);
//		for (Node n2 : g.nodes) {
//			//System.out.println(((Col)n2.getAttr("Col")).getColor());
//			if(sn.id == n2.id)
//				continue;
//			Edge e = g.addEdge(sn.id, n2.id);
//			e.putAttr("Type", 0);
//			e.putAttr("Value", betaSn);
//		}
		mg.addGraph(g);
		return mg;
	}

	public static MultiGraph SIS(double beta, double delta, int N) {
		MultiGraph mg = new MultiGraph();
		Graph gb = new Graph(false);
		gb.putAttr("type", "beta");
		Edge eb = gb.addEdge(0, 1);
		eb.putAttr("value", beta);
		mg.addGraph(gb);
		
		Graph gd = new Graph(false);
		gd.putAttr("type", "delta");
		Edge ed = gb.addEdge(1, 0);
		ed.putAttr("value", delta);
		mg.addGraph(gd);
		
		mg.getAddNode(0).color = new Color(1f,1f,1f);
		mg.getAddNode(1).color = new Color(0f,0f,0f);
		
		return mg;
	}
	
	public static MultiGraph MultipleDiseaseCompetitionGraph(int M) {
		double[] beta = {0.001, 0.01};
		String[] types = {"ss", "ba", "sg"};
		MultiGraph mg = new MultiGraph();
		Graph g = new Graph(false);
		g.putAttr("type", "beta");
		Node ns = new Node(0);
		g.nodes.add(ns);
		for (int i = 1; i <= M; i++) {
			Node n = new Node(i);
			g.nodes.add(n);
			g.addEdge(0, i).putAttr("value", beta[i-1]);
		}
		mg.addGraph(g);
		for(MultiNode mn: mg.nodes)
			mn.putAttr("type", types[mn.id]);
		return mg;
	}

	public static MultiGraph SIIS(double[] beta, double[] delta) {
		MultiGraph mg = new MultiGraph();
		Graph gb = new Graph(false);
		gb.putAttr("type", "beta");
		Edge eb = gb.addEdge(0, 1);
		eb.putAttr("value", beta[0]);
		Edge eb2 = gb.addEdge(0, 2);
		eb2.putAttr("value", beta[1]);
		mg.addGraph(gb);
		
		Graph gd = new Graph(false);
		gd.putAttr("type", "delta");
		Edge ed = gd.addEdge(1, 0);
		ed.putAttr("value", delta[0]);
		Edge ed2 = gd.addEdge(2, 0);
		ed2.putAttr("value", delta[1]);
		mg.addGraph(gd);
		
		mg.getAddNode(0).color = new Color(1f,1f,1f);
		mg.getAddNode(1).color = new Color(0f,0f,0f);
		mg.getAddNode(2).color = new Color(.5f,.5f,.5f);
		
		return mg;
	}
	
	public static MultiGraph I123(double[] beta){
		MultiGraph mg = new MultiGraph();
		Graph gb = new Graph(false);
		gb.putAttr("type", "beta");
		Edge eb = gb.addEdge(0, 1);
		eb.putAttr("value", beta[0]);
		Edge eb2 = gb.addEdge(1, 2);
		eb2.putAttr("value", beta[1]);
		Edge eb3 = gb.addEdge(2, 0);
		eb3.putAttr("value", beta[2]);
		mg.addGraph(gb);
				
		return mg;
	}

	public static MultiGraph SIISbeta(double[] beta) {
		MultiGraph mg = new MultiGraph();
		Graph gb = new Graph(false);
		gb.putAttr("type", "beta");
		Edge eb = gb.addEdge(0, 1);
		eb.putAttr("value", beta[0]);
		Edge eb2 = gb.addEdge(0, 2);
		eb2.putAttr("value", beta[1]);
		Edge eb3 = gb.addEdge(1, 0);
		eb3.putAttr("value", beta[2]);
		Edge eb4 = gb.addEdge(2, 0);
		eb4.putAttr("value", beta[3]);
		mg.addGraph(gb);
				
		return mg;
	}
	
}
