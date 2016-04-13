package cs.contagion;

import java.awt.Container;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.security.SecureRandom;
import java.util.Date;
import java.util.Random;

import javax.imageio.ImageIO;
import javax.swing.JFrame;

import org.math.plot.Plot2DPanel;

import cs.graph.Edge;
import cs.graph.Graph;
import cs.graph.MultiGraph;
import cs.graph.MultiNode;
import cs.graph.Node;
import cs.graph.datasets.Files;

public class Utils {
	

	public static void SaveDiseaseGraph(boolean saveMg, MultiGraph mg, String folderName, String sampleName) throws IOException {
		if(!saveMg) return;
		File f = new File("results/data/" + folderName + "/" + sampleName + "mg.txt");
		if (!f.exists()) {
			f.getParentFile().mkdirs();
			f.createNewFile();
		}
		BufferedWriter bw = new BufferedWriter(new FileWriter(f));
		bw.write("Number of states: " + mg.nodes.size() + "\n");
		bw.write("Beta matrix:\n");
		for (int i = 0; i < mg.nodes.size(); i++) {
			int ind = 0;
			for (Edge e : mg.nodes.get(i).outEdges("beta")) {
				while (e.dest.id > ind) {
					bw.write(0 + ", ");
					ind++;
				}
				bw.write((Double) e.getAttr("value") + (e.dest.id == (mg.nodes.size() - 1) ? ";\n" : ", "));
				ind++;
			}
			while (mg.nodes.size() > ind) {
				ind++;
				bw.write(0 + (ind == mg.nodes.size() ? ";\n" : ", "));
			}
		}
		bw.write("Delta matrix:\n");
		for (int i = 0; i < mg.nodes.size(); i++) {
			int ind = 0;
			for (Edge e : mg.nodes.get(i).outEdges("delta")) {
				while (e.dest.id > ind) {
					bw.write(0 + ", ");
					ind++;
				}
				bw.write((Double) e.getAttr("value") + (e.dest.id == (mg.nodes.size() - 1) ? ";\n" : ", "));
				ind++;
			}
			while (mg.nodes.size() > ind) {
				ind++;
				bw.write(0 + (ind == mg.nodes.size() ? ";\n" : ", "));
			}
		}
		bw.close();
	}

	public static void AssignNodeDiseases(int assignment, Graph g, MultiGraph mg, float[][] colors, boolean saveMg, String filename) throws IOException {
		int M = mg.nodes.size();
		int N = g.nodes.size();
		int L = (int) Math.sqrt(N);
		Random ran = new SecureRandom();
		if (assignment == 0) { // ???
			for (Node n : g.nodes) {
				if (ran.nextDouble() <= 0.9995) {
					n.putAttr("curCol", mg.nodes.get(mg.nodes.size() - 1));
					continue;
				}
				colors[n.id] = ((MultiNode) n.getAttr("curCol")).getColor();
			}
		} else if (assignment == 1) { // assigned to closest source (language)
			for (Node n : g.nodes) {
				int minInd = -1;
				double minDist = Double.MAX_VALUE;
				for (int i = 0; i < M; i++) {
					int[] coords = (int[]) mg.nodes.get(i).getAttr("sourceCoords");
					double dist = Math.sqrt(Math.pow(Math.abs(n.id % L - coords[0]), 2) + Math.pow(Math.abs(n.id / L - coords[1]), 2));
					if (dist < minDist) {
						minDist = dist;
						minInd = i;
					}
				}
				n.putAttr("curCol", mg.nodes.get(minInd));
				colors[n.id] = ((MultiNode) n.getAttr("curCol")).getColor();
			}
		} else if (assignment == 2) { // only sources assigned (language with
			// susceptibility)
			for (Node n : g.nodes) {
				n.putAttr("curCol", mg.nodes.get(M));
				for (int i = 0; i < M; i++) {
					int[] coords = (int[]) mg.nodes.get(i).getAttr("sourceCoords");
					if (n.id == coords[0] % L + coords[1] * L) {
						n.putAttr("curCol", mg.nodes.get(i));
						break;
					}
				}
				colors[n.id] = ((MultiNode) n.getAttr("curCol")).getColor();
			}
		} else if (assignment == 3) // random assignment
			for (Node n : g.nodes) {
				n.putAttr("curCol", mg.nodes.get(ran.nextInt(mg.nodes.size())));
				colors[n.id] = ((MultiNode) n.getAttr("curCol")).getColor();
			}
		else if (assignment == 4) {
			// int ch = ran.nextInt(N);
			// int[] ds = {0,1,2,0};
			// n.putAttr("curCol", gCols.nodes.get(ds[n.id]));
			// n.putAttr("curCol", nColS);
			// if (n.id == ch)
			// n.putAttr("curCol", gCols.nodes.get(0));
			// n.putAttr("curCol", col[((n.id % L) + (n.id / L)) % 2]);//
			// grid-like
		} else if (assignment == 5) { // from file - 3 diseases
			if (filename != "") {
				BufferedReader br = new BufferedReader(new FileReader(filename));
				String line = "";
				while ((line = br.readLine()) != null) {
					String[] strs = line.split(" ");
					if (g.nodes.size() <= Integer.parseInt(strs[0]))
						break;
					g.nodes.get(Integer.parseInt(strs[0])).putAttr("curCol", mg.nodes.get(Integer.parseInt(strs[1])));
				}
				br.close();
			}
		}

		if (saveMg && (filename != "") && (assignment == 3)) {
			File f = new File(filename);
			if (!f.exists()) {
				f.getParentFile().mkdirs();
				f.createNewFile();
			}
			BufferedWriter bw = new BufferedWriter(new FileWriter(f));
			for (Node n : g.nodes) {
				bw.write(n.id + " " + ((MultiNode) n.getAttr("curCol")).id);
				bw.newLine();
			}
			bw.close();
		}
	}

	public static MultiGraph SetDiseaseGraph(int setDiseaseGraph, int M, int N, double[] beta, double[] delta) {
		MultiGraph mg = null;
		if (setDiseaseGraph == 0)
			mg = MultiGraph.FullMeshDiseaseToyGraph(M, N);
		else if (setDiseaseGraph == 1) {
			double beta_mean = 0.15;
			mg = MultiGraph.createFullMeshDiseaseGraph(M, N, beta_mean, 0.0001);
		} else if (setDiseaseGraph == 2) {
			mg = MultiGraph.LanguageSpread(N);
		} else if (setDiseaseGraph == 3) {
			mg = MultiGraph.SIS(0.5, 0.1, N);
		} else if (setDiseaseGraph == 4) {
			mg = MultiGraph.MultipleDiseaseCompetitionGraph(M);
		} else if (setDiseaseGraph == 5) {
			mg = MultiGraph.SIIS(beta, delta);
		} else if (setDiseaseGraph == 6) {
			mg = MultiGraph.I123(beta);
		} else if (setDiseaseGraph == 7) {
			mg = MultiGraph.SIISbeta(beta);
		}
		return mg;
	}

	public static Graph SetGraph(int setGraph, int N, String[] args) throws IOException {
		int L = (int) Math.sqrt(N);
		Graph g = null;
		String filename = args[0] + "/" + args[1] + "graph.txt";
		if (setGraph == 0)
			g = Graph.readGraphEdgeList(filename, true);
		else if (setGraph == 1)
			g = Graph.readGraphGML(Files.files[17], true, false);
		else if (setGraph == 2)
			g = Graph.createSquareGridGraph(L, 3, true);// 0 - random 2-hop, 1 -
		// small world, 2 -
		// random node, 3 - grid
		// g.saveGraphEdgeList(filename);
		// System.out.println(g.nedges * 1.0 / N);
		// g.saveGraphEdgeList("results/data/"+sampleName+"graph.txt");
		N = g.nodes.size();
		return g;
	}

	public static void saveMatlab(boolean saveMatlab, MultiGraph mg, Graph g, int N, String folderName, String sampleName, int iter) throws IOException{
		if (!saveMatlab) return;
		double[] ci = new double[mg.nodes.size()];
		int[][] st = new int[N][mg.nodes.size()];
		for (Node n : g.nodes) {
			ci[((MultiNode) n.getAttr("nextCol")).id] += 1.0 / N;
			st[n.id][((MultiNode) n.getAttr("nextCol")).id] = 1;
		}
		BufferedWriter w = new BufferedWriter(new FileWriter("results/data/" + folderName + "/" + sampleName + iter + "-sim-mv.txt", true));
		String ww = "";
		for (int i = 0; i < mg.nodes.size(); i++) {
			if (i > 0)
				ww += ", ";
			ww += ci[i] + "";
		}
		ww += "\n";

		w.write(ww);
		w.close();
	}
	
	public static void plotJava(boolean plotJava, MultiGraph mg, Graph g, int it, int N, double[] countInfected, Plot2DPanel plot){
		int stepPlotJava = 50;
		if (plotJava && it % stepPlotJava == 0) {
			if (it < stepPlotJava) {
				for (Node n : g.nodes) {
					countInfected[((Node) n.getAttr("nextCol")).id] += 1.0 / N;
					// countInfected[((Col) ((Node)
					// n.getAttr("nextCol")).getAttr("Col")).id] += 1.0 / N;
				}
				return;
			}
			double[] countInfectedNew = new double[mg.nodes.size()];
			for (Node n : g.nodes) {
				countInfectedNew[((Node) n.getAttr("nextCol")).id] += 1.0 / N;
				// countInfectedNew[((Col) ((Node)
				// n.getAttr("nextCol")).getAttr("Col")).id] += 1.0 / N;
			}
			double[] xx = new double[] { it - stepPlotJava, it };
			for (int i = 0; i < mg.nodes.size(); i++) {
				double[] yy = new double[] { countInfected[i], countInfectedNew[i] };
				plot.addLinePlot(String.valueOf(i), (mg.nodes.get(i)).color, new double[] { xx[0], yy[0] }, new double[] { xx[1], yy[1] });
				// plot.addLinePlot(String.valueOf(i),((Col) ((Node)
				// gCols.nodes.get(i)).getAttr("Col")).getColor(), new
				// double[] { xx[0], yy[0] },new double[] { xx[1], yy[1]
				// });// yy);//Arrays.copyOfRange(countInfected[i],0,it+1));
			}
			countInfected = countInfectedNew.clone();
			plot.setFixedBounds(1, 0, 1);
			plot.repaint();
		}
	}

	public static void saveImages(boolean saveImages, int it, Container contentPane, String folderName, String sampleName, int iter) throws IOException{
		int stepSaveImages = 10;
		if (saveImages && it % stepSaveImages == 0 && it > 0) {
			BufferedImage img = new BufferedImage(contentPane.getWidth(), contentPane.getHeight(), BufferedImage.TYPE_INT_RGB);
			contentPane.paint(img.getGraphics());
			File f = new File("results/images/" + folderName + "/" + sampleName + iter + "-" + it / stepSaveImages + ".png");
			if (!f.exists()) {
				f.getParentFile().mkdirs();
				f.createNewFile();
			}
			ImageIO.write(img, "png", f);
		}
	}

	public static Plot2DPanel setJavaPlot(boolean plotJava){
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
		return plot;
	}
	
	public static String Toc(Date tic) {
		Date toc = new Date();
		long diff = toc.getTime() - tic.getTime();
		int mins = (int) Math.floor((diff / 1000.0) / 60.0);
		int secs = (int) Math.floor(diff / 1000.0 - mins * 60.0);
		return ((mins < 10) ? ("0" + mins) : mins) + ":" + ((secs < 10) ? ("0" + secs) : secs);
	}
	
}
