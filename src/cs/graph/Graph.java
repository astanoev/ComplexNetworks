package cs.graph;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.RandomAccessFile;
import java.io.StringReader;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

public class Graph {

	public boolean symmetric;

	public Graph(boolean symmetric) {
		this.symmetric = symmetric;
	}

	public Graph(Graph other) {
		this.symmetric = other.symmetric;
		this.nedges = other.nedges;
		this.nodes = other.nodes;
		this.attr = other.attr;
	}

	public ArrayList<Node> nodes = new ArrayList<Node>();

	public int nedges = 0;

	public static final boolean PRINT_PROGRESS = true;

	private HashMap<String, Object> attr = new HashMap<String, Object>();

	public Node addNode(int nodeID) {
		int index = Collections.binarySearch(nodes, nodeID);
		if (index >= 0) {
			System.out.println("The node with ID " + nodeID + " already exists");
			return null;
		} else {
			Node tmp = new Node(nodeID);
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

	public List<Node> commonNodes(List<Edge> a, List<Edge> b) {
		List<Node> commonNodes = new ArrayList<Node>();
		int i = 0, j = 0;
		while (i < a.size() && j < b.size()) {
			if (a.get(i).dest.id == b.get(j).dest.id) {
				commonNodes.add(a.get(i).dest);
				i++;
				j++;
			} else if (a.get(i).dest.id > b.get(j).dest.id) {
				j++;
			} else {
				i++;
			}
		}
		return commonNodes;
	}

	public Edge addEdge(int srcID, int destID) {
		int indexSrc = Collections.binarySearch(nodes, srcID);
		Node srcNode = null;
		if (indexSrc < 0) {
			srcNode = addNode(srcID);
		} else {
			srcNode = nodes.get(indexSrc);
		}
		Node destNode = null;
		int indexDest = Collections.binarySearch(nodes, destID);
		if (indexDest < 0) {
			destNode = addNode(destID);
		} else {
			destNode = nodes.get(indexDest);
		}
		Edge edge = srcNode.addEdge(destNode);
		if (edge != null) {
			destNode.addInEdge(edge);
			nedges++;
		}

		if (symmetric && edge != null) {
			Edge edge2 = addEdge(destID, srcID);
			// if(edge.attr!=null)
			// edge2.attr = edge.attr;
		}
		return edge;
	}

	public boolean hasEdge(int srcID, int destID) {
		int indexSrc = Collections.binarySearch(nodes, srcID);
		Node srcNode = null;
		if (indexSrc < 0) {
			return false;
		} else {
			srcNode = nodes.get(indexSrc);
		}
		Node destNode = new Node(destID);
		return srcNode.hasEdge(destNode);
	}

	public Edge getEdge(int srcID, int destID) {
		int indexSrc = Collections.binarySearch(nodes, srcID);
		Node srcNode = null;
		if (indexSrc < 0) {
			return null;
		} else {
			srcNode = nodes.get(indexSrc);
		}
		Node destNode = null;
		int indexDest = Collections.binarySearch(nodes, destID);
		if (indexDest < 0) {
			return null;
		} else {
			destNode = nodes.get(indexDest);
		}
		return srcNode.getEdge(destNode);
	}

	public void saveGraphEdgeList(String filename) throws IOException {
		System.out.println(filename);
		File f = new File(filename);
		if (!f.exists()) {
			f.getParentFile().mkdirs();
			f.createNewFile();
		}
		BufferedWriter w = new BufferedWriter(new FileWriter(filename));
		for (Node node : nodes) {
			for (Edge edge : node.outEdges) {
				w.write(node.id + " " + edge.dest.id + "\n");
			}
		}
		w.close();
	}

	public static Graph readGraphPajek(String filename, boolean makeSymmetric, boolean printConsole) throws NumberFormatException, IOException {

		InputStream is = null;
		try {
			// url = new URL(getCodeBase(), filename);
			is = Graph.class.getClassLoader().getResourceAsStream(filename);
		} catch (Exception e) {
			if (printConsole)
				System.out.println(filename + " : Malformed URL ");
			return null;
		}
		BufferedReader r = null;
		try {
			// r = new BufferedReader(new InputStreamReader(url.openStream()));
			if (filename.startsWith("0"))
				r = new BufferedReader(new StringReader(filename.substring(1)));
			else
				r = new BufferedReader(new InputStreamReader(is));
		} catch (Exception ex) {
			if (printConsole)
				System.out.println("rabbit " + ex.getMessage());
			return null;
		}
		String line = null;
		Graph g = new Graph(makeSymmetric);

		// int lineno = 0;
		// double oneProcEdges = -1;
		// double oneProcVertices = -1;
		int state = 0;
		// int prevProc = 0;
		// read all the edges
		while ((line = r.readLine()) != null) {

			line = line.trim();
			String[] ss = line.split("[ \t]+");

			if (ss[0].startsWith("*Vertices")) {
				state = 1;
				// oneProcVertices = Integer.parseInt(ss[1]) / 100.0;
				continue;
			} else if (line.startsWith("*Arcs")) {
				state = 2;
				continue;
			} else if (line.startsWith("*Edges")) {
				state = 3;
				continue;
			}
			if (state == 1) {

				int id = Integer.parseInt(ss[0]);
				Node node = g.addNode(id);
				node.putAttr("label", ss[1].replace('"', '\0'));
				try {
					node.putAttr("x", ss[2]);
					node.putAttr("y", ss[3]);
					node.putAttr("z", ss[4]);
				} catch (Exception ex) {
					if (printConsole)
						System.out.println("No coordinates in dataset");
				}
			} else if (state == 3) {
				int a = Integer.parseInt(ss[0]);
				int b = Integer.parseInt(ss[1]);
				Edge edge = g.addEdge(a, b);
				double value = 1.0;
				if (ss.length >= 3) { // weigths...
					value = Double.parseDouble(ss[2]);
				}
				if (edge == null) {
					if (printConsole)
						System.out.println("rabbit");
					continue;
				}
				edge.putAttr("value", value);
				if (g.symmetric) {
					Edge edgeR = g.addEdge(b, a);
					edgeR.putAttr("value", value);
				}

			}
		}
		if (printConsole)
			System.out.println("nodes/edges: " + g.nodes.size() + "/" + g.nedges);
		return g;

	}

	public static Graph readGraphEdgeList(String filename, boolean makeSymmetric) throws NumberFormatException, IOException {
		return readGraphEdgeList2(filename, makeSymmetric, false);
	}

	public static Graph readGraphEdgeList2(String filename, boolean makeSymmetric, boolean print) throws IOException {
		Graph g = new Graph(makeSymmetric);
		File f = new File(filename);
		if (!f.exists()) {
			f.getParentFile().mkdirs();
			f.createNewFile();
		}
		FileChannel inChannel = new RandomAccessFile(filename, "r").getChannel();
		MappedByteBuffer buffer = inChannel.map(FileChannel.MapMode.READ_ONLY, 0, inChannel.size());
		buffer.load();
		String line = "";
		for (int i = 0; i < buffer.limit(); i++) {
			char ch = (char) buffer.get();
			if (ch != '\n')
				line = line + ch;
			else {
				if (line.startsWith("#"))
					continue;
				line.trim();
				String[] ss = line.trim().split(" +");
				line = "";
				Edge edge = null;
				double value = 1.0;
				int a = Integer.parseInt(ss[0]);
				int b = Integer.parseInt(ss[1]);
				edge = g.addEdge(a, b);
				if (ss.length >= 3) {// weigths...
					value = Double.parseDouble(ss[2]);
				}
				if (edge == null) {
					if (print)
						System.out.println("rabbit");
					continue;
				}
				edge.putAttr("value", value);
				if (g.symmetric) {
					Edge edgeR = g.addEdge(b, a);
					edgeR.putAttr("value", value);
				}
			}
		}
		buffer.clear();
		inChannel.close();
		return g;
	}

	public static Graph readGraphEdgeList(String filename, boolean makeSymmetric, boolean print) throws NumberFormatException, IOException {

		InputStream is = null;
		try {
			// url = new URL(getCodeBase(), filename);
			// is = new Graph(true).getClass().getResourceAsStream(filename);
			is = new FileInputStream(filename);
		} catch (Exception e) {
			if (print)
				System.out.println(filename + " : Malformed URL ");
			return null;
		}
		BufferedReader r = null;
		try {
			// r = new BufferedReader(new InputStreamReader(url.openStream()));
			if (filename.startsWith("0"))
				r = new BufferedReader(new StringReader(filename.substring(1)));
			else
				r = new BufferedReader(new InputStreamReader(is));
		} catch (Exception ex) {
			if (print)
				System.out.println("rabbit " + ex.getMessage());
			return null;
		}
		Graph g = new Graph(makeSymmetric);
		int lineno = 0;
		double oneProcEdges = -1;
		int prevProc = 0;
		String line = "";
		// read all the edges
		while ((line = r.readLine()) != null) {
			if (line.startsWith("#"))
				continue;
			line = line.trim();
			String[] ss = line.split("[ \t]+");
			if (ss.length == 1) {
				oneProcEdges = Integer.parseInt(ss[0]) / 100.0;
				continue;
			}
			lineno++;
			if (PRINT_PROGRESS && oneProcEdges != -1 && print) {
				int tmp = (int) Math.floor(lineno / oneProcEdges);
				if (tmp > prevProc) {
					prevProc = tmp;
					System.out.println(prevProc + "%");
				}
			}
			int a = 0;
			int b = 0;
			Edge edge = null;
			double value = 1.0;
			// try{
			a = Integer.parseInt(ss[0]);
			b = Integer.parseInt(ss[1]);
			edge = g.addEdge(a, b);
			if (ss.length >= 3) {// weigths...
				value = Double.parseDouble(ss[2]);
			}
			// }
			// catch(Exception ex){
			// return null;
			// }
			if (edge == null) {
				if (print)
					System.out.println("rabbit");
				continue;
			}
			edge.putAttr("value", value);
			if (g.symmetric) {
				Edge edgeR = g.addEdge(b, a);
				edgeR.putAttr("value", value);
			}

		}
		if (print)
			System.out.println("nodes/edges: " + g.nodes.size() + "/" + g.nedges);
		return g;
	}

	public static Graph readGraphGML(String filename, boolean makeSymmetric, boolean printConsole) throws IOException {
		InputStream is = null;
		try {
			// url = new URL(getCodeBase(), filename);
			// is = Graph.class.getClassLoader().getResourceAsStream(filename);
			is = new FileInputStream("datasets/gml/" + filename);
		} catch (Exception e) {
			if (printConsole)
				System.out.println(filename + " : Malformed URL ");
			return null;
		}
		BufferedReader r = null;
		try {
			// r = new BufferedReader(new InputStreamReader(url.openStream()));
			if (filename.startsWith("0"))
				r = new BufferedReader(new StringReader(filename.substring(1)));
			else
				r = new BufferedReader(new InputStreamReader(is));
		} catch (Exception ex) {
			if (printConsole)
				System.out.println("rabbit " + ex.getMessage());
			return null;
		}
		String line = null;
		Graph g = new Graph(makeSymmetric);

		int id = -1;
		String label = null;
		double value = Double.MAX_VALUE;
		int source = 0;
		int target = 0;
		// whereami tells in which section we are in; node=0;edge=1
		int whereami = -1;
		while ((line = r.readLine()) != null) {
			line = line.trim();
			if (line.startsWith("]")) {
				if (whereami == 0) {
					whereami = -1;
					Node node = g.addNode(id);
					if (label != null) {
						node.putAttr("label", label);
						label = null;
					}
					if (value != Double.MAX_VALUE) {
						node.putAttr("value", value);
						value = Double.MAX_VALUE;
					}
				} else if (whereami == 1) {
					whereami = -1;
					Edge edge = g.addEdge(source, target);

					if (value == Double.MAX_VALUE)
						value = 1.0;
					if (edge == null) {
						if (printConsole)
							System.out.println("double link!!!");
						continue;
					}
					edge.putAttr("value", value);
					if (g.symmetric == true) {
						Edge edgeR = g.addEdge(target, source);
						edgeR.putAttr("value", value);
					}
					value = Double.MAX_VALUE;
				}
			} else if (line.startsWith("node")) {
				whereami = 0;
			} else if (line.startsWith("edge")) {
				whereami = 1;
			} else if (line.startsWith("label")) {
				label = line.split("[ \t]+")[1].replace("\"", "");
			} else if (line.startsWith("id")) {
				id = Integer.parseInt(line.split("[ \t]+")[1]);
			} else if (line.startsWith("source")) {
				source = Integer.parseInt(line.split("[ \t]+")[1]);
			} else if (line.startsWith("target")) {
				target = Integer.parseInt(line.split("[ \t]+")[1]);
			} else if (line.startsWith("value")) {
				value = Double.parseDouble(line.split("[ \t]+")[1]);
			}
		}
		if (printConsole)
			System.out.println("nodes/edges: " + g.nodes.size() + "/" + g.nedges);
		return g;

	}

	public static void makeZeroIndexFile(String filename) throws IOException {
		File f = new File(filename);
		File fTemp = new File(filename + "-temp");
		BufferedReader br = new BufferedReader(new FileReader(f));
		BufferedWriter bw = new BufferedWriter(new FileWriter(fTemp));
		String line = "";
		while ((line = br.readLine()) != null) {
			String[] ss = line.trim().split(" +");
			int a = Integer.parseInt(ss[0]) - 1;
			int b = Integer.parseInt(ss[1]) - 1;
			String outLine = a + " " + b;
			for (int i = 2; i < ss.length; i++)
				outLine += " " + ss[i];
			bw.write(outLine + "\n");
		}
		br.close();
		bw.close();
		f.delete();
		fTemp.renameTo(f);
	}

	public static Graph createSquareGridGraph(int width, int smallWorld, boolean Torus) {
		int N = width * width;
		double alfa_a = 2;
		double[] p = new double[width - 1];
		double[] pp = new double[width];
		int[] cc = new int[width - 1];
		int[] ccc = new int[width - 1];

		if (smallWorld == 1) {
			pp[0] = 0;
			for (int i = 2; i <= width; i++) {
				if (i == width / 2) {
					cc[i - 2] = 2 * (width - 1);
				} else if (i < width / 2) {
					cc[i - 2] = 4 * i;
				} else if (i == width) {
					cc[i - 2] = 1;
				} else {
					cc[i - 2] = 4 * (width - i);
				}
				p[i - 2] = cc[i - 2] * Math.pow(i, -alfa_a);
				pp[i - 1] = pp[i - 2] + p[i - 2];
			}
		}

		Random ran = new SecureRandom();
		Graph g = new Graph(true);
		for (int i = 0; i < N; i++) {
			if ((i + 1) % width != 0) {
				g.addEdge(i, i + 1);
				g.addEdge(i + 1, i);
			} else if (Torus) {
				g.addEdge(i, i + 1 - width);
				g.addEdge(i + 1 - width, i);
			}
			if (i / width + 1 != width) {
				g.addEdge(i, i + width);
				g.addEdge(i + width, i);
			} else if (Torus) {
				g.addEdge(i, i % width);
				g.addEdge(i % width, i);
			}
			if ((smallWorld == 1) && (ran.nextDouble() <= 0.05)) {// small
																	// world,
																	// alpha = 2
				double ranG = ran.nextDouble() * pp[width - 2];
				int ii = binarySearch(pp, ranG, 1, width - 1) - 1;
				int jj = ran.nextInt(cc[ii]);
				int x = (int) (4.0 * (ii + 2) / cc[ii] * jj - (ii + 2));
				if (jj > cc[ii] / 2.0)
					x = (int) (-4.0 * (ii + 2) / cc[ii] * jj + 3.0 * (ii + 2));
				int y = (int) (-4.0 * (ii + 2) / cc[ii] * jj + 2.0 * (ii + 2));
				if (jj <= cc[ii] / 4.0)
					y = (int) (4.0 * (ii + 2) / cc[ii] * jj);
				else if (jj > 3 * cc[ii] / 4.0)
					y = (int) (4.0 * (ii + 2) / cc[ii] * jj - 4.0 * (ii + 2));
				ccc[ii]++;
				g.addEdge(i, (i + x + y * width + N) % N);
				// g.addEdge((i + x + y * L + N) % N,i);
			} else if (smallWorld == 2) {// choose additional random node
				int neig = i;
				while (true)
					if ((neig = ran.nextInt(N)) != i) {
						g.addEdge(i, neig);
						break;
					}
			} else if (smallWorld == 0) {// choose additional random 2-hop
											// neighbor
				int jj = ran.nextInt(8);
				int ii = 0;
				int x = (int) (4.0 * (ii + 2) / 8 * jj - (ii + 2));
				if (jj > 8 / 2.0)
					x = (int) (-4.0 * (ii + 2) / 8 * jj + 3.0 * (ii + 2));
				int y = (int) (-4.0 * (ii + 2) / 8 * jj + 2.0 * (ii + 2));
				if (jj <= 8 / 4.0)
					y = (int) (4.0 * (ii + 2) / 8 * jj);
				else if (jj > 3 * 8 / 4.0)
					y = (int) (4.0 * (ii + 2) / 8 * jj - 4.0 * (ii + 2));
				g.addEdge(i, (i + x + y * width + N) % N);
			}
		}
		return g;
	}

	private static int binarySearch(double[] pp, double ranG, int l, int r) {
		int i = (r + l) / 2;
		if ((ranG <= pp[i]) && (ranG > pp[i - 1]))
			return i;
		else if (ranG > pp[i])
			return binarySearch(pp, ranG, i + 1, r);
		else
			return binarySearch(pp, ranG, l, i - 1);
	}

	public static Graph createBarabasiAlbertNetwork(int N, int m0, int m) {
		Graph g = new Graph(true);
		int totalWeights = m0 * (m0 - 1);
		int[] nodeWeights = new int[N];
		Random rnd = new SecureRandom();
		for (int i = 0; i < m0 - 1; i++) {
			g.addNode(i);
			for (int j = i + 1; j < m0; j++) {
				if (rnd.nextDouble() <= m * 1.0 / (m0 - 1)) {
					g.addEdge(i, j).putAttr("value", 1.0);
					//g.addEdge(j, i).putAttr("value", 1.0);
				}
			}
			nodeWeights[i] = g.nodes.get(i).outEdges.size();
		}
		g.addNode(m0 - 1);
		nodeWeights[m0 - 1] = g.nodes.get(m0 - 1).outEdges.size();
		for (int i = m0; i < N; i++) {
			int[] cn = chooseNeighbors(i, nodeWeights, m, totalWeights);
			for (int j = 0; j < m; j++) {
				nodeWeights[cn[j]]++;
				g.addEdge(i, cn[j]).putAttr("value", 1.0);
				//g.addEdge(cn[j], i).putAttr("value", 1.0);
			}
			nodeWeights[i] = m;
			totalWeights += 2 * m;
		}
		int[] nodeDegrees = new int[N];
		for (Node n : g.nodes)
			nodeDegrees[n.outEdges.size()]++;
		return g;
	}

	private static int[] chooseNeighbors(int i, int[] nodeWeights, int m, int totalWeights) {
		int[] cn = new int[m];
		boolean[] cb = new boolean[i];
		Random rnd = new SecureRandom();
		for (int j = 0; j < m; j++) {
			double sum = 0;
			for (int k = 0; k < i; k++) {
				if (cb[k])
					continue;
				double ch = rnd.nextDouble();
				sum += nodeWeights[k] * 1.0 / totalWeights;
				if (ch < sum) {
					cn[j] = k;
					cb[k] = true;
					totalWeights -= nodeWeights[k];
					break;
				}
			}
		}
		return cn;
	}

	public static Graph createWSN(int n, double range) {
		Graph g = new Graph(true);
		Random rnd = new SecureRandom();
		while (true) {
			g = new Graph(true);
			for (int i = 0; i < n; i++) {
				Node node = g.addNode(i);
				node.putAttr("locX", rnd.nextDouble());
				node.putAttr("locY", rnd.nextDouble());
				for (int j = 0; j < i; j++) {
					Node node2 = g.nodes.get(j);
					if (distance(node, node2) <= range)
						g.addEdge(node.id, node2.id);
				}
			}
			ArrayList<Integer> nodes = new ArrayList<Integer>();
			LinkedList<Integer> q = new LinkedList<Integer>();
			q.push(0);
			while (!q.isEmpty()) {
				int id = q.pop();
				nodes.add(id);
				Node node = g.nodes.get(id);
				for (Edge e : node.outEdges)
					if (!nodes.contains(e.dest.id) && !q.contains(e.dest.id))
						q.push(e.dest.id);
			}
			boolean ok = true;
			for (Node node : g.nodes)
				// if (node.outEdges.size() == 0) {
				if (!nodes.contains(node.id)) {
					ok = false;
					break;
				}
			if (ok)
				break;
		}
		return g;
	}

	private static double distance(Node node1, Node node2) {
		return Math.sqrt(((Double) node2.getAttr("locX") - (Double) node1.getAttr("locX")) * ((Double) node2.getAttr("locX") - (Double) node1.getAttr("locX"))
				+ ((Double) node2.getAttr("locY") - (Double) node1.getAttr("locY")) * ((Double) node2.getAttr("locY") - (Double) node1.getAttr("locY")));
	}

	public void measurements1Hop(int nSamples, double vd) {
		Random rnd = new SecureRandom();
		for (Node node : this.nodes) {
			for (Edge e : node.outEdges) {
				double dist = distance(node, e.dest);
				double sigma = vd * dist;
				double[] mmts = new double[nSamples];
				double mean = 0;
				double var = 0;
				for (int i = 0; i < nSamples; i++) {
					mmts[i] = dist + sigma * rnd.nextGaussian();
					mean += mmts[i];
				}
				mean /= nSamples;
				for (int i = 0; i < nSamples; i++) {
					var += (mean - mmts[i])*(mean - mmts[i]);
				}
				var /= (nSamples-1.0);
				e.putAttr("mean", mean);
				e.putAttr("var", var);
			}
		}
	}

	public Graph measurements2Hop(double r) {
		Graph g2 = new Graph(true);
		for (Node node : this.nodes) {
			g2.addNode(node.id);
			for (Edge e1 : node.outEdges) {
				double mean1 = (Double) e1.getAttr("mean");
				double var1 = (Double) e1.getAttr("var");
				for (Edge e2 : e1.dest.outEdges) {
					if (e2.dest.id == node.id)
						continue;
					if (this.getEdge(node.id, e2.dest.id) != null)
						continue;
					double mean2 = (Double) e2.getAttr("mean");
					double var2 = (Double) e2.getAttr("var");
					// calculate with 2-hop approximation - motion update step
					double cosAlfa = -Math.sqrt((r * r - (mean1 - mean2) * (mean1 - mean2)) / (4 * mean1 * mean2));
					double mean = Math.sqrt(mean1 * mean1 + mean2 * mean2 - 2 * mean1 * mean2 * cosAlfa);
					double var = var1 + var2;

					Edge e = g2.getEdge(node.id, e2.dest.id);
					if (e != null) {
						if (e.getAttr("mean") != null) {
							double meanOld = (Double) e.getAttr("mean");
							double varOld = (Double) e.getAttr("var");
							// combine with previous measurements - measurement
							// update step
							double meanNew = (meanOld / varOld + mean / var) / (1 / varOld + 1 / var);
							double varNew = 1.0 / (1 / varOld + 1 / var);
							e.putAttr("mean", meanNew);
							e.putAttr("var", varNew);
						} else {
							e.putAttr("mean", mean);
							e.putAttr("var", var);
						}
					} else {
						e = g2.addEdge(node.id, e2.dest.id);
						e.putAttr("mean", mean);
						e.putAttr("var", var);
					}
				}
			}
		}
		return g2;
	}
}
