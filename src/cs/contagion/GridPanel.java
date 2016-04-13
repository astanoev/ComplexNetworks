package cs.contagion;

import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;
import javax.swing.JFrame;
import javax.swing.JPanel;

public class GridPanel extends JPanel {
	private static final long serialVersionUID = 6102003947415373483L;
	private int L;
	public float[][] colors;
	private int pixelWidth;

	public GridPanel(int pixelWidth, int L, float[][] colors) {
		this.pixelWidth = pixelWidth;
		this.L = L;
		this.colors = colors;
	}

	public void paintComponent(Graphics g) {
		super.paintComponent(g);

		for (int i = 0; i < L * L; i++) {
			g.setColor(new Color(colors[i][0], colors[i][1], colors[i][2]));
			g.fillRect((i % L) * pixelWidth, (i / L) * pixelWidth, pixelWidth, pixelWidth);
		}
	}
	
	public void setDrawGraph(boolean drawGraph, boolean saveImages, JFrame frame, Container contentPane, int L, int width, String folderName, String sampleName, int iter) throws IOException{
		if (!drawGraph) return;
		frame.setTitle("GridRect");
		frame.setLocation(100, 100);
		// frame.setSize(L * width, L * width);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		Dimension d = new Dimension(L * width, L * width);
		contentPane.setPreferredSize(d);
		frame.pack();
		frame.setVisible(true);
		contentPane.add(this);
		contentPane.repaint();
		if (saveImages) {
			BufferedImage img = new BufferedImage(contentPane.getWidth(), contentPane.getHeight(), BufferedImage.TYPE_INT_RGB);
			contentPane.paint(img.getGraphics());
			File f = new File("results/images/" + folderName + "/" + sampleName + iter + "-0.png");
			if (!f.exists()) {
				f.getParentFile().mkdirs();
				f.createNewFile();
			}
			ImageIO.write(img, "png", f);
		}
	}
}
