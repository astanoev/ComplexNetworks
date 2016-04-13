package cs.contagion;

import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.image.BufferedImage;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.imageio.ImageIO;
import javax.swing.JFrame;

import cs.contagion.GridPanel;
import cs.graph.Edge;
import cs.graph.Graph;
import cs.graph.MultiGraph;
import cs.graph.MultiNode;
import cs.graph.Node;

import org.math.plot.Plot2DPanel;

public class MultiplexSpread {

	public static void main(String[] args) throws IOException, InterruptedException {
		// parameters
		boolean plotJava = false;
		boolean saveMatlab = false;
		boolean saveImages = false;
		String sampleName = "LanguageSpread3.2.1-";
		int M = 2;
		int L = 100; // power of 2, because of the parallel jobs
		int N = L * L;
		int width = 4;
		float[][] colors = new float[N][3];
		int setGraph = 0;
		int setDiseaseGraph = 4;
		int nodeAssignment = 3;

		MultiGraph g = SetGraph(setGraph, N);
		MultiGraph mg = Utils.SetDiseaseGraph(setDiseaseGraph, M, N, null, null);
		AssignNodeDiseases(nodeAssignment, g, mg, colors);

		JFrame frame = new JFrame();
		frame.setTitle("GridRect");
		frame.setLocation(100, 100);

		// frame.setSize(L * width, L * width);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		Container contentPane = frame.getContentPane();
		Dimension d = new Dimension(L * width, L * width);
		contentPane.setPreferredSize(d);
		frame.pack();
		frame.setVisible(true);
		GridPanel gp = new GridPanel(width, L, colors);
		contentPane.add(gp);

		Plot2DPanel plot = new Plot2DPanel();
		JFrame framePlot = new JFrame("Number of infected nodes");
		if (plotJava) {
			framePlot.setContentPane(plot);
			framePlot.setVisible(true);
			framePlot.setSize(1000, 600);
			framePlot.setLocation(750, 100);
			plot.setAxisLabel(0, "Time step");
			plot.setAxisLabel(1, "Percentage of infected nodes");
			plot.setFixedBounds(1, 0, 1);
			plot.repaint();
		}
		double[] countInfected = new double[mg.nodes.size()];

		int nrOfProcessors = Runtime.getRuntime().availableProcessors();
		int tasksPerProcessor = 2;// L/nrOfProcessors;
		//ExecutorService eservice = Executors.newFixedThreadPool(nrOfProcessors * tasksPerProcessor);
		ExecutorService eservice = Executors.newFixedThreadPool(N);
		CompletionService<List<MultiNode>> cservice = new ExecutorCompletionService<List<MultiNode>>(
				eservice);

		BufferedWriter w = null;
		if(saveMatlab){
			w = new BufferedWriter(new FileWriter("results/data/"
					+ sampleName + "ci.txt"));
			w.close();
		}
		int maxIt = 5000;
		int step = 2;//Math.max(L / (nrOfProcessors * tasksPerProcessor), 1);
		long begTest = 0;
		for (int it = 0; it < maxIt; it++) {
			frame.setTitle("GridRect  " + it);

			for (int index = 0; index < L * L; index += L * step)
				cservice.submit(new GridTask(g, index, index + L * step, mg));

			// List<Node> GridTaskResult;
			// g.nodes = new ArrayList<Node>();
			for (int index = 0; index < L * L; index += L * step) {
				try {
					cservice.take().get();
					// GridTaskResult = cservice.take().get();
					// g.nodes.addAll(GridTaskResult);
				} catch (InterruptedException e) {
				} catch (ExecutionException e) {
				}
			}

//			 //for debugging the parallel job
//			 GridTask gt = new GridTask(g, 0, N, mg);
//			 gt.call();

			for (Node n : g.nodes) {
				n.putAttr("curCol", (MultiNode) n.getAttr("nextCol"));
				colors[n.id] = ((MultiNode) n.getAttr("nextCol")).getColor();
			}
			// Collections.sort((List)g.nodes);
			gp.colors = colors;
			contentPane.repaint();
			//Thread.sleep(500);
			if (it % 1000 == 0) {
				Double secs = new Double(
						(new java.util.Date().getTime() - begTest) * 0.001);
				System.out.println("run time " + secs + " secs");
				begTest = new java.util.Date().getTime();
			}

			// for plot in matlab
			if (saveMatlab) {
				double[] ci = new double[mg.nodes.size()];
				int[][] st = new int[N][mg.nodes.size()];
				for (Node n : g.nodes) {
					ci[((Node) n.getAttr("nextCol")).id] += 1.0 / N;
					// ci[((Col) ((Node)
					// n.getAttr("nextCol")).getAttr("Col")).id] += 1.0 / N;
					st[n.id][((Node) n.getAttr("nextCol")).id] = 1;
					// st[n.id][((Col) ((Node)
					// n.getAttr("nextCol")).getAttr("Col")).id] = 1;
				}
				w = new BufferedWriter(new FileWriter("results/data/"
						+ sampleName + "ci.txt", true));
				for (int i = 0; i < mg.nodes.size(); i++) {
					if (i > 0)
						w.write(", ");
					w.write(ci[i] + "");
				}
				w.newLine();
				w.close();
				BufferedWriter w2 = new BufferedWriter(new FileWriter(
						"results/data/" + sampleName + it + "-st.txt", true));
				for (int i = 0; i < N; i++) {
					for (int j = 0; j < mg.nodes.size(); j++) {
						if (j > 0)
							w2.write(", ");
						w2.write(st[i][j] + "");
					}
					w2.newLine();
				}
				w2.close();
			}

			int stepPlotJava = 50;
			if (plotJava && it % stepPlotJava == 0) {
				if (it < stepPlotJava) {
					for (Node n : g.nodes) {
						countInfected[((Node) n.getAttr("nextCol")).id] += 1.0 / N;
						// countInfected[((Col) ((Node)
						// n.getAttr("nextCol")).getAttr("Col")).id] += 1.0 / N;
					}
					continue;
				}
				double[] countInfectedNew = new double[mg.nodes.size()];
				for (Node n : g.nodes) {
					countInfectedNew[((Node) n.getAttr("nextCol")).id] += 1.0 / N;
					// countInfectedNew[((Col) ((Node)
					// n.getAttr("nextCol")).getAttr("Col")).id] += 1.0 / N;
				}
				double[] xx = new double[] { it - stepPlotJava, it };
				for (int i = 0; i < mg.nodes.size(); i++) {
					double[] yy = new double[] { countInfected[i],
							countInfectedNew[i] };
					plot.addLinePlot(String.valueOf(i),
							FloatArrayToColor((float[]) (mg.nodes.get(i))
									.getAttr("color")), new double[] { xx[0],
									yy[0] }, new double[] { xx[1], yy[1] });
					// plot.addLinePlot(String.valueOf(i),((Col) ((Node)
					// gCols.nodes.get(i)).getAttr("Col")).getColor(), new
					// double[] { xx[0], yy[0] },new double[] { xx[1], yy[1]
					// });// yy);//Arrays.copyOfRange(countInfected[i],0,it+1));
				}
				countInfected = countInfectedNew.clone();
				plot.setFixedBounds(1, 0, 1);
				plot.repaint();
			}
			int stepSaveImages = 10;
			if (saveImages && it % stepSaveImages == 0 && it > 0) {
				BufferedImage img = new BufferedImage(contentPane.getWidth(),
						contentPane.getHeight(), BufferedImage.TYPE_INT_RGB);
				contentPane.paint(img.getGraphics());
				ImageIO.write(img, "png", new File("results/images/"
						+ sampleName + it / stepSaveImages + ".png"));
			}
		}
	}

	private static void AssignNodeDiseases(int nodeAssignment, MultiGraph g,
			MultiGraph mg, float[][] colors) {
		Random rn = new SecureRandom();
		for(MultiNode n: g.nodes){
			n.putAttr("curCol", mg.nodes.get(0));
			colors[n.id] = mg.nodes.get(0).getColor();
		}
		for(int i=0;i<1;i++){
		MultiNode mn1 = g.nodes.get(rn.nextInt(g.nodes.size()));
		mn1.putAttr("curCol", mg.nodes.get(1));
		colors[mn1.id] = mg.nodes.get(1).getColor();
		MultiNode mn2 = g.nodes.get(rn.nextInt(g.nodes.size()));
		mn2.putAttr("curCol", mg.nodes.get(2));
		colors[mn2.id] = mg.nodes.get(2).getColor();
		}
	}

	private static MultiGraph SetGraph(int setGraph, int N) throws IOException {
		MultiGraph mg = new MultiGraph();
		if(setGraph==0){
			Graph g1 = Graph.createBarabasiAlbertNetwork(N, 10, 2);
			System.out.println(g1.nedges);
			g1.putAttr("type", "ba");
			//g1.saveGraphEdgeList("datasets/edgeformat/ba-sample-003.txt");
			Graph g2 = Graph.createSquareGridGraph((int)Math.sqrt(N), 1, false);
			System.out.println(g2.nedges);
			g2.putAttr("type", "sg");
			mg.addGraph(g1);
			mg.addGraph(g2);
		}	
		return mg;
	}

	private static Color FloatArrayToColor(float[] attr) {
		return new Color(attr[0], attr[1], attr[2]);
	}

	public static class GridTask implements Callable<List<MultiNode>>{
		private MultiGraph g;
		private int i_start;
		private int i_end;
		private MultiGraph mg;

		public GridTask() {
		}

		public GridTask(MultiGraph g2, int i_start, int i_end, MultiGraph mg) {
			this.g = g2;
			this.i_start = i_start;
			this.i_end = i_end;
			this.mg = mg;
		}

		public List<MultiNode> call() {

			Random ran = new SecureRandom();
			for (MultiNode n : g.nodes.subList(i_start, i_end)) {
				try{
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
					String[] dis = {"ba", "sg"};
					for(String ds : dis)
						for (Edge e : n.getNode(ds).outEdges) {
							if((String)((MultiNode)((MultiNode)e.dest.getAttr("superNode")).getAttr("curCol")).getAttr("type")!=ds)
								continue;
							Edge ed = mn.getNode("beta").getEdge(((MultiNode)((MultiNode)e.dest.getAttr("superNode")).getAttr("curCol")).getNode("beta"));
							if (ed != null)
								if (ran.nextDouble() <= (Double)ed.getAttr("value"))
									infectiveNeighbors.add(e.dest);
						}
					if (infectiveNeighbors.size() > 0) 
						n.putAttr("nextCol",(MultiNode) ((MultiNode)infectiveNeighbors.get(ran.nextInt(infectiveNeighbors.size())).getAttr("superNode")).getAttr("curCol"));
				}
				catch(Exception ex){
					System.out.println("error");
				}
			}

			return g.nodes.subList(i_start, i_end);
		}
	}

}
