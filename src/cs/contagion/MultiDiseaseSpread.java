package cs.contagion;

import java.awt.Container;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import javax.swing.JFrame;

import cs.contagion.GridPanel;
import cs.graph.Edge;
import cs.graph.Graph;
import cs.graph.MultiGraph;
import cs.graph.MultiNode;
import cs.graph.Node;
import cs.contagion.Utils;

import org.math.plot.Plot2DPanel;

public class MultiDiseaseSpread {

	public static void main(String[] args) throws InterruptedException, NumberFormatException, IOException {
		mainFunction(0);
	}
	
	public static void batchSimulations() throws InterruptedException, NumberFormatException, IOException {
		// for executing multiple simulations and saving them and so..
		int initIter = 0;
		int maxIter = 500;
		ThreadPoolExecutor executor = new ThreadPoolExecutor(7, maxIter - initIter, Long.MAX_VALUE, TimeUnit.MINUTES, new ArrayBlockingQueue<Runnable>(maxIter - initIter));
		for (int i = initIter; i < maxIter; i++) {
			final int j = i;
			executor.execute(new Runnable() {
				@Override
				public void run() {
					Date tic = new Date();
					System.out.print("Iteration " + j + " underway...");
					try {
						// i = mainFunction(i);
						mainFunction(j);
					} catch (Exception e) {
						e.printStackTrace();
					}
					System.out.println(" finished in " + Utils.Toc(tic) + " minutes");
				}
			});
		}
		executor.shutdown();
		executor.awaitTermination(Long.MAX_VALUE, TimeUnit.MINUTES);
		System.exit(0);
	}

	public static int mainFunction(int iter) throws InterruptedException, NumberFormatException, IOException {
		// parameters
		int M = 3; // number of states (contagions)
		int L = 256; // grid width (power of 2, when parallelizing)
		int N = L * L; // number of nodes (==number of grid tiles)
		float[][] colors = new float[N][3];
		// params for saving files..
		String sampleName = "power-";
		String folderName = "Graph17";
		String[] args = { "datasets/edgeformat/" + folderName + "/", sampleName };// iter+"-"};
		// params for setting the network, state transition network and initial values
		int setGraph = 2;
		int setDiseaseGraph = 0;
		int nodeAssignment = 3;
		boolean saveMg = false;

		// PlosPaperColorMap(N, M, setGraph, setDiseaseGraph, nodeAssignment, folderName, sampleName, args);

		// create complex network
		Graph g = Utils.SetGraph(setGraph, N, args);
		N = g.nodes.size();
		
		// create network of state(contagion) transitions 
		double[] delta = new double[] { 0.075, 0.05 };
		double[] beta = new double[] { 0.09, 0.05, 0.15 };
		MultiGraph mg = Utils.SetDiseaseGraph(setDiseaseGraph, M, N, beta, delta);
		M = mg.nodes.size();
		Utils.SaveDiseaseGraph(saveMg, mg, folderName, sampleName);
		
		// assign the nodes of the network with initial values (states/contagions)
		Utils.AssignNodeDiseases(nodeAssignment, g, mg, colors, saveMg, args[0] + args[1] + "ic.txt");

		Simulation(N, L, colors, g, mg, folderName, sampleName, iter);
		
		// calcDeterministicApproximation(g, mg, folderName, sampleName, true);
		
		// g = SetGraph(setGraph, N, args);
		// N = g.nodes.size();
		// mg = SetDiseaseGraph(setDiseaseGraph, M, N, beta, delta);
		// M = mg.nodes.size();
		// colors = new float[N][3];
		// AssignNodeDiseases(nodeAssignment, g, mg, colors, args[0] + args[1] +
		// "ic.txt");
		//
		// calcMonteCarloApproximation(g, mg, folderName, sampleName, iter);
		// calcDeterministicReal(g, mg, folderName, sampleName, true);

		return iter;
	}

	public static void Simulation(int N, int L, float[][] colors, Graph g, MultiGraph mg, String folderName, String sampleName, int iter) throws IOException {

		boolean plotJava = false;
		boolean saveMatlab = false;
		boolean saveImages = false;
		boolean drawGraph = true;
		int pixelWidth = 1;

		for (Node n : g.nodes)
			colors[n.id] = ((MultiNode) n.getAttr("curCol")).getColor();

		JFrame frame = new JFrame();
		Container contentPane = frame.getContentPane();
		GridPanel gp = new GridPanel(pixelWidth, L, colors);
		gp.setDrawGraph(drawGraph, saveImages, frame, contentPane, L, pixelWidth, folderName, sampleName, iter);
		Plot2DPanel plot = Utils.setJavaPlot(plotJava);
		
		double[] countInfected = new double[mg.nodes.size()];

		int nrOfProcessors = Runtime.getRuntime().availableProcessors();
		int tasksPerProcessor = 4;// L/nrOfProcessors;
		ExecutorService eservice = Executors.newFixedThreadPool(N); //nrOfProcessors * tasksPerProcessor
		CompletionService<List<Node>> cservice = new ExecutorCompletionService<List<Node>>(eservice);

		BufferedWriter w = null;
		if (saveMatlab) {
			File f = new File("results/data/" + folderName + "/" + sampleName + iter + "-sim-mv.txt");
			if (!f.exists()) {
				f.getParentFile().mkdirs();
				f.createNewFile();
			}
			w = new BufferedWriter(new FileWriter(f));
			w.close();
		}
		int maxIt = 1000;
		int step = Math.max(N / (nrOfProcessors * tasksPerProcessor), 1);
		for (int it = 0; it < maxIt; it++) {
			frame.setTitle("GridRect  " + it);

			 for (int index = 0; index < N; index += step)
			 cservice.submit(new GridTask(g, index, index + Math.min(step,N-index), mg));
			
			 // List<Node> GridTaskResult;
			 // g.nodes = new ArrayList<Node>();
			 for (int index = 0; index < N; index += step) {
			 try {
				 cservice.take().get();
				 // GridTaskResult = cservice.take().get();
				 // g.nodes.addAll(GridTaskResult);
			 } catch (InterruptedException e) {
			 } catch (ExecutionException e) {
			 }
			 }

			// for debugging the parallel job
//			GridTask gt = new GridTask(g, 0, N, mg);
//			gt.call();

			for (Node n : g.nodes) {
				n.putAttr("curCol", (MultiNode) n.getAttr("nextCol"));
				colors[n.id] = ((MultiNode) n.getAttr("nextCol")).getColor();
			}
			if (drawGraph) {
				gp.colors = colors;
				contentPane.repaint();
			}
			
			Utils.saveMatlab(saveMatlab, mg, g, N, folderName, sampleName, iter);
			Utils.plotJava(plotJava, mg, g, it, N, countInfected, plot);
			Utils.saveImages(saveImages, it, contentPane, folderName, sampleName, iter);
		}
	}

	public static void PlosPaperColorMap(int N, int M, int setGraph, int setDiseaseGraph, int nodeAssignment, String folderName, String sampleName, String[] args)
			throws IOException, InterruptedException {
		// for beta1-beta2 colormap of the paper
		boolean saveMg = false;
		float[][] colors = new float[N][3];
		double[] delta = new double[] { 0.075, 0.05 };
		double[] beta = new double[4];
		beta[2] = 0.075;
		beta[3] = 0.05;
		double stepBeta = 0.01;
		for (beta[0] = 0; beta[0] <= 1; beta[0] = Math.round(100 * (beta[0] + stepBeta)) / 100.0)
			for (beta[1] = 0; beta[1] <= 1; beta[1] = Math.round(100 * (beta[1] + stepBeta)) / 100.0) {

				Graph g = Utils.SetGraph(setGraph, N, args);
				MultiGraph mg = Utils.SetDiseaseGraph(setDiseaseGraph, M, N, beta, delta);
				Utils.AssignNodeDiseases(nodeAssignment, g, mg, colors, saveMg, args[0] + args[1] + "ic.txt");
				N = g.nodes.size();
				M = mg.nodes.size();
				System.out.println("BETA 1: " + beta[0] + " , BETA 2: " + beta[1]);
				calcDeterministicApproximation(g, mg, folderName, sampleName + beta[0] + "-" + beta[1] + "-", false);
				g = Utils.SetGraph(setGraph, N, args);
				mg = Utils.SetDiseaseGraph(setDiseaseGraph, M, N, beta, delta);
				colors = new float[N][3];
				Utils.AssignNodeDiseases(nodeAssignment, g, mg, colors, saveMg, args[0] + args[1] + "ic.txt");
				N = g.nodes.size();
				M = mg.nodes.size();
				calcDeterministicReal(g, mg, folderName, sampleName + beta[0] + "-" + beta[1] + "-", false);
				// calcMonteCarloApproximation(g, mg, folderName,
				// sampleName, iter);
			}
	}

	public static void PlosPaperColorMap2(int N, int M, final int setGraph, final int setDiseaseGraph, final int nodeAssignment, final String folderName, final String sampleName,
			final String[] args) throws IOException, InterruptedException {
		// for beta1-beta2 colormap of the paper
		boolean saveMg = false;
		float[][] colors = new float[N][3];
		double[] delta = new double[] { 0.075, 0.05 };
		double[] beta = new double[4];
		beta[2] = 0.075;
		beta[3] = 0.05;
		double stepBeta = 0.01;
		int nrComb = 2 * ((int) Math.round(1.0 / stepBeta + 1)) * ((int) Math.round(1.0 / stepBeta + 1));
		ThreadPoolExecutor executor = new ThreadPoolExecutor(1, nrComb, Long.MAX_VALUE, TimeUnit.MINUTES, new ArrayBlockingQueue<Runnable>(nrComb));
		for (beta[0] = 0; beta[0] <= 1; beta[0] = Math.round(100 * (beta[0] + stepBeta)) / 100.0)
			for (beta[1] = 0; beta[1] <= 1; beta[1] = Math.round(100 * (beta[1] + stepBeta)) / 100.0) {
				final double[] betaN = beta;
				final double[] deltaN = delta;
				System.out.println("BETA 1: " + betaN[0] + " , BETA 2: " + betaN[1]);

				final Graph g = Utils.SetGraph(setGraph, N, args);
				// System.out.println("BETA 1: " + betaN[0]
				// + " , BETA 2: " + betaN[1]);
				final MultiGraph mg = Utils.SetDiseaseGraph(setDiseaseGraph, M, N, betaN, deltaN);
				Utils.AssignNodeDiseases(nodeAssignment, g, mg, colors, saveMg, args[0] + args[1] + "ic.txt");
				N = g.nodes.size();
				M = mg.nodes.size();
				// System.out.println("BETA 1: " + betaN[0]
				// + " , BETA 2: " + betaN[1]);
				executor.execute(new Runnable() {
					@Override
					public void run() {
						try {
							calcDeterministicApproximation(g, mg, folderName, sampleName + betaN[0] + "-" + betaN[1] + "-", false);
						} catch (Exception ex) {
							ex.printStackTrace();
							System.out.println(ex.getMessage());
						}
					}
				});
				final Graph g2 = Utils.SetGraph(setGraph, N, args);
				final MultiGraph mg2 = Utils.SetDiseaseGraph(setDiseaseGraph, M, N, betaN, deltaN);
				colors = new float[N][3];
				Utils.AssignNodeDiseases(nodeAssignment, g2, mg2, colors, saveMg, args[0] + args[1] + "ic.txt");
				N = g2.nodes.size();
				M = mg2.nodes.size();
				executor.execute(new Runnable() {
					@Override
					public void run() {
						try {
							calcDeterministicReal(g2, mg2, folderName, sampleName + betaN[0] + "-" + betaN[1] + "-", false);
							// calcMonteCarloApproximation(g, mg, folderName,
							// sampleName, iter);
						} catch (Exception ex) {
							ex.printStackTrace();
							System.out.println(ex.getMessage());
						}
					}
				});
			}
		executor.shutdown();
		executor.awaitTermination(Long.MAX_VALUE, TimeUnit.MINUTES);
		System.exit(0);
	}

	private static void calcDeterministicApproximation(Graph g, MultiGraph mg, String folderName, String sampleName, boolean printIterations) throws IOException {
		File f = new File("results/data/" + folderName + "/" + sampleName + "appr-mean.txt");
		File f2 = new File("results/data/" + folderName + "/" + sampleName + "appr-std.txt");
		if (f.exists())
			f.delete();
		f.getParentFile().mkdirs();
		f.createNewFile();
		if (f2.exists())
			f2.delete();
		f2.getParentFile().mkdirs();
		f2.createNewFile();
		double[][] p = new double[g.nodes.size()][mg.nodes.size()];
		for (Node n : g.nodes)
			p[n.id][((MultiNode) n.getAttr("curCol")).id] = 1.0;
		int maxit = 50000;
		int checkIt = 500;
		double[] pCheck = new double[mg.nodes.size()];
		for (int it = 1; it <= maxit; it++) {
			if (printIterations)
				System.out.print("Iteration " + it + " underway... ");
			double[][] pNew = new double[g.nodes.size()][mg.nodes.size()];
			double[][] pB = new double[g.nodes.size()][mg.nodes.size()];
			for (Node n : g.nodes)
				for (MultiNode mnK : mg.nodes)
					for (Edge e : mnK.outEdges("beta"))
						pB[n.id][mnK.id] += ((Double) e.getAttr("value")) * p[n.id][e.dest.id];
			for (Node n : g.nodes) {
				for (MultiNode mnK : mg.nodes) {
					double deltaSum = 0.0;
					for (Edge e : mnK.outEdges("delta")) {
						pNew[n.id][e.dest.id] += p[n.id][mnK.id] * ((Double) e.getAttr("value"));
						deltaSum += ((Double) e.getAttr("value"));
					}
					double den = 0.0;
					double gK0 = 0.0;
					for (Edge e : n.outEdges) {
						den += pB[e.dest.id][mnK.id] / (1.0 - 0.5 * pB[e.dest.id][mnK.id]);
						if (1.0 > pB[e.dest.id][mnK.id])
							gK0 += Math.log(1.0 - pB[e.dest.id][mnK.id]);
						else
							gK0 += Math.log(0.0);
					}
					gK0 = Math.exp(gK0);
					for (Edge e : mnK.outEdges("beta")) {
						double num = 0.0;
						for (Edge e2 : n.outEdges) {
							num += p[e2.dest.id][e.dest.id] / (1.0 - 0.5 * pB[e2.dest.id][mnK.id]);
						}
						if (den == 0)
							break;
						// System.out.println("err");
						pNew[n.id][e.dest.id] += p[n.id][mnK.id] * (1 - deltaSum) * ((Double) e.getAttr("value")) * num / den * (1 - gK0);
					}
					pNew[n.id][mnK.id] += p[n.id][mnK.id] * (1 - deltaSum) * gK0;
				}
			}
			for (Node n : g.nodes)
				for (MultiNode mn : mg.nodes)
					if (pNew[n.id][mn.id] <= 1e-10)
						pNew[n.id][mn.id] = 0.0;
			for (Node n : g.nodes) {
				double sum = 0.0;
				for (MultiNode mnK : mg.nodes)
					sum += pNew[n.id][mnK.id];
				for (MultiNode mnK : mg.nodes)
					pNew[n.id][mnK.id] /= sum;
			}
			for (Node n : g.nodes)
				p[n.id] = Arrays.copyOf(pNew[n.id], pNew[n.id].length);
			double[] pp = new double[mg.nodes.size()];
			BufferedWriter bw = new BufferedWriter(new FileWriter(f, true));
			for (MultiNode mn : mg.nodes) {
				for (Node n : g.nodes)
					pp[mn.id] += p[n.id][mn.id];
				pp[mn.id] = pp[mn.id] / g.nodes.size();
				if (printIterations)
					bw.write(pp[mn.id] + ((mn.id == (mg.nodes.size() - 1)) ? "\n" : ", "));
			}
			bw.close();
			double[] pp2 = new double[mg.nodes.size()];
			BufferedWriter bw2 = new BufferedWriter(new FileWriter(f2, true));
			for (MultiNode mn : mg.nodes) {
				for (Node n : g.nodes)
					pp2[mn.id] += Math.abs(p[n.id][mn.id]-pp[mn.id]);
				pp2[mn.id] = pp2[mn.id] / g.nodes.size();
				if (printIterations)
					bw2.write(pp2[mn.id] + ((mn.id == (mg.nodes.size() - 1)) ? "\n" : ", "));
			}
			bw2.close();
			if (printIterations)
				System.out.println("done!");
			if (it % checkIt == 0) {
				double norm = 0.0;
				for (MultiNode mn : mg.nodes)
					norm += (pp[mn.id] - pCheck[mn.id]) * (pp[mn.id] - pCheck[mn.id]);
				norm = Math.sqrt(norm);
				if (!printIterations && ((norm <= 1E-20) || (it == maxit))) {
					bw = new BufferedWriter(new FileWriter(f, true));
					for (MultiNode mn : mg.nodes)
						bw.write(pp[mn.id] + ((mn.id == (mg.nodes.size() - 1)) ? "\n" : ", "));
					bw.close();
					break;
				}
				pCheck = Arrays.copyOf(pp, pp.length);
			}
		}
	}

	private static void calcDeterministicReal(Graph g, MultiGraph mg, String folderName, String sampleName, boolean printIterations) throws IOException {
		File f = new File("results/data/" + folderName + "/" + sampleName + "real.txt");
		if (f.exists())
			f.delete();
		f.getParentFile().mkdirs();
		f.createNewFile();
		double[][] p = new double[g.nodes.size()][mg.nodes.size()];
		for (Node n : g.nodes)
			p[n.id][((MultiNode) n.getAttr("curCol")).id] = 1.0;
		int maxit = 50000;
		int checkIt = 500;
		double[] pCheck = new double[mg.nodes.size()];
		for (int it = 1; it <= maxit; it++) {
			if (printIterations)
				System.out.print("Iteration " + it + " underway... ");
			double[][] pNew = new double[g.nodes.size()][mg.nodes.size()];
			double[][] pB = new double[g.nodes.size()][mg.nodes.size()];
			for (Node n : g.nodes)
				for (MultiNode mnK : mg.nodes)
					for (Edge e : mnK.outEdges("beta"))
						pB[n.id][mnK.id] += ((Double) e.getAttr("value")) * p[n.id][e.dest.id];
			for (Node n : g.nodes) {
				for (MultiNode mnK : mg.nodes) {
					double deltaSum = 0.0;
					for (Edge e : mnK.outEdges("delta")) {
						pNew[n.id][e.dest.id] += p[n.id][mnK.id] * ((Double) e.getAttr("value"));
						deltaSum += ((Double) e.getAttr("value"));
					}
					double gK0 = 0.0;
					for (Edge e : n.outEdges)
						if (1.0 > pB[e.dest.id][mnK.id])
							gK0 += Math.log(1.0 - pB[e.dest.id][mnK.id]);
						else
							gK0 += Math.log(0.0);
					gK0 = Math.exp(gK0);
					double[] fk = fK(n, mnK, p, pB);
					int nBe = 0;
					for (Edge e : mnK.outEdges("beta")) {
						pNew[n.id][e.dest.id] += p[n.id][mnK.id] * (1 - deltaSum) * fk[nBe++];
					}
					pNew[n.id][mnK.id] += p[n.id][mnK.id] * (1 - deltaSum) * gK0;
				}
			}
			for (Node n : g.nodes)
				for (MultiNode mn : mg.nodes)
					if (pNew[n.id][mn.id] <= 1e-10)
						pNew[n.id][mn.id] = 0.0;
			for (Node n : g.nodes) {
				double sum = 0.0;
				for (MultiNode mnK : mg.nodes)
					sum += pNew[n.id][mnK.id];
				for (MultiNode mnK : mg.nodes)
					pNew[n.id][mnK.id] /= sum;
			}
			for (Node n : g.nodes)
				p[n.id] = Arrays.copyOf(pNew[n.id], pNew[n.id].length);
			double[] pp = new double[mg.nodes.size()];
			BufferedWriter bw = new BufferedWriter(new FileWriter(f, true));
			for (MultiNode mn : mg.nodes) {
				for (Node n : g.nodes)
					pp[mn.id] += p[n.id][mn.id];
				pp[mn.id] = pp[mn.id] / g.nodes.size();
				if (printIterations)
					bw.write(pp[mn.id] + ((mn.id == (mg.nodes.size() - 1)) ? "\n" : ", "));
			}
			bw.close();
			if (printIterations)
				System.out.println("done!");
			if (it % checkIt == 0) {
				double norm = 0.0;
				for (MultiNode mn : mg.nodes)
					norm += (pp[mn.id] - pCheck[mn.id]) * (pp[mn.id] - pCheck[mn.id]);
				norm = Math.sqrt(norm);
				if (!printIterations && ((norm <= 1E-20) || (it == maxit))) {
					bw = new BufferedWriter(new FileWriter(f, true));
					for (MultiNode mn : mg.nodes)
						bw.write(pp[mn.id] + ((mn.id == (mg.nodes.size() - 1)) ? "\n" : ", "));
					bw.close();
					break;
				}
				pCheck = Arrays.copyOf(pp, pp.length);
			}
		}
	}

	private static double[] fK(Node n, MultiNode mnK, double[][] p, double[][] pB) {
		double[] f = new double[mnK.outEdges("beta").size()];
		for (int i = 1; i < Math.pow(2, n.outEdges.size()); i++) { // each event
			String event = Integer.toBinaryString(i);
			while (event.length() < n.outEdges.size())
				event = "0" + event;
			ArrayList<Node> infNeighbors = new ArrayList<Node>(); // infective
			// neighbors
			for (int j = 0; j < event.length(); j++)
				if (event.charAt(j) == '1')
					infNeighbors.add(n.outEdges.get(j).dest);
			// double[] gK0 = new double[mnK.outEdges("beta").size()]; //gK0 for
			// the non-infective neighbors
			double gK0 = 0.0; // gK0 for the non-infective neighbors
			for (int j = 0; j < n.outEdges.size(); j++) {
				if (event.charAt(j) == '0')
					// for(int l=0; l<mnK.outEdges("beta").size();l++)
					// gK0[l] += Math.log(1.0 -
					// pB[n.outEdges.get(j).dest.id][mnK.outEdges("beta").get(l).dest.id]);
					if (1.0 > pB[n.outEdges.get(j).dest.id][mnK.id])
						gK0 += Math.log(1.0 - pB[n.outEdges.get(j).dest.id][mnK.id]);
					else
						gK0 += Math.log(0.0);
			}

			for (int k = 0; k < Math.pow(mnK.outEdges("beta").size(), infNeighbors.size()); k++) { // each
																									// infective
																									// configuration
				String conf = Integer.toString(k, mnK.outEdges("beta").size());
				while (conf.length() < infNeighbors.size())
					conf = "0" + conf;
				double gK1 = 0.0;
				for (int j = 0; j < infNeighbors.size(); j++) {
					Edge e = mnK.outEdges("beta").get(Character.getNumericValue(conf.charAt(j)));
					gK1 += Math.log(p[infNeighbors.get(j).id][e.dest.id]) + Math.log((Double) e.getAttr("value"));
					if (gK1 == Double.NEGATIVE_INFINITY)
						break;
				}
				for (int j = 0; j < infNeighbors.size(); j++)
					f[Character.getNumericValue(conf.charAt(j))] += 1.0 / infNeighbors.size() * Math.exp(gK1 + gK0);
				// f[Character.getNumericValue(conf.charAt(j))] +=
				// 1.0/infNeighbors.size()*Math.exp(gK1)*Math.exp(gK0[Character.getNumericValue(conf.charAt(j))]);
			}
		}
		return f;
	}

	@SuppressWarnings("unused")
	private static void calcMonteCarloApproximation(Graph g, MultiGraph mg, String folderName, String sampleName, int iter) throws IOException {
		File f = new File("results/data/" + folderName + "/" + sampleName + iter + "-mc.txt");
		if (!f.exists()) {
			f.getParentFile().mkdirs();
			f.createNewFile();
		} else
			f.delete();
		double[][] p = new double[g.nodes.size()][mg.nodes.size()];
		for (Node n : g.nodes)
			p[n.id][((MultiNode) n.getAttr("curCol")).id] = 1.0;
		Random rnd = new SecureRandom();
		int maxit = 1000;
		for (int it = 1; it <= maxit; it++) {
			// System.out.print("Iteration "+it+" underway... ");
			double[][] pNew = new double[g.nodes.size()][mg.nodes.size()];
			double[][] pB = new double[g.nodes.size()][mg.nodes.size()];
			for (Node n : g.nodes)
				for (MultiNode mnK : mg.nodes)
					for (Edge e : mnK.outEdges("beta"))
						pB[n.id][mnK.id] += ((Double) e.getAttr("value")) * p[n.id][e.dest.id];
			for (Node n : g.nodes) {
				for (MultiNode mnK : mg.nodes) {
					double deltaSum = 0.0;
					for (Edge e : mnK.outEdges("delta")) {
						pNew[n.id][e.dest.id] += p[n.id][mnK.id] * ((Double) e.getAttr("value"));
						deltaSum += ((Double) e.getAttr("value"));
					}
					double den = 0.0;
					double gK0 = 0.0;
					for (Edge e : n.outEdges) {
						den += pB[e.dest.id][mnK.id] / (1.0 - 0.5 * pB[e.dest.id][mnK.id]);
						if (1.0 > pB[e.dest.id][mnK.id])
							gK0 += Math.log(1.0 - pB[e.dest.id][mnK.id]);
						else
							gK0 += Math.log(0.0);
					}
					gK0 = Math.exp(gK0);
					for (Edge e : mnK.outEdges("beta")) {
						double num = 0.0;
						for (Edge e2 : n.outEdges) {
							num += p[e2.dest.id][e.dest.id] / (1.0 - 0.5 * pB[e2.dest.id][mnK.id]);
						}
						if (den == 0)
							break;
						// System.out.println("err");
						pNew[n.id][e.dest.id] += p[n.id][mnK.id] * (1 - deltaSum) * ((Double) e.getAttr("value")) * num / den * (1 - gK0);
					}
					pNew[n.id][mnK.id] += p[n.id][mnK.id] * (1 - deltaSum) * gK0;
				}
			}
			for (Node n : g.nodes)
				for (MultiNode mn : mg.nodes)
					if (pNew[n.id][mn.id] <= 1e-10)
						pNew[n.id][mn.id] = 0.0;
			for (Node n : g.nodes) {
				double sum = 0.0;
				for (MultiNode mnK : mg.nodes)
					sum += pNew[n.id][mnK.id];
				for (MultiNode mnK : mg.nodes)
					pNew[n.id][mnK.id] /= sum;
			}
			p = new double[g.nodes.size()][mg.nodes.size()];
			for (Node n : g.nodes) {
				double prob = rnd.nextDouble();
				double cum = 0.0;
				for (int i = 0; i < mg.nodes.size(); i++) {
					if (prob <= cum + pNew[n.id][mg.nodes.get(i).id]) {
						p[n.id][mg.nodes.get(i).id] = 1.0;
						break;
					}
					cum += pNew[n.id][mg.nodes.get(i).id];
				}
			}
			BufferedWriter bw = new BufferedWriter(new FileWriter(f, true));
			for (MultiNode mn : mg.nodes) {
				double pp = 0.0;
				for (Node n : g.nodes)
					pp += p[n.id][mn.id];
				bw.write(pp / g.nodes.size() + ((mn.id == (mg.nodes.size() - 1)) ? "\n" : ", "));
			}
			bw.close();
		}
	}

}
