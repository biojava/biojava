/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */


package org.biojava.stats.svm.tools;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Cursor;
import java.awt.FlowLayout;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Paint;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.geom.AffineTransform;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.Collections;
import java.util.Iterator;
import java.util.Set;

import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JRadioButton;

import org.biojava.stats.svm.PolynomialKernel;
import org.biojava.stats.svm.RadialBaseKernel;
import org.biojava.stats.svm.SMOTrainer;
import org.biojava.stats.svm.SVMClassifierModel;
import org.biojava.stats.svm.SVMKernel;
import org.biojava.stats.svm.SVMTarget;
import org.biojava.stats.svm.SimpleSVMTarget;

/**
 * A simple toy example that allows you to put points on a canvas, and find a
 * polynomial hyperplane to seperate them.
 *
 * @author Ewan Birney
 * @author Matthew Pocock
 * @author Thomas Down
 * @author Michael L Heurer
 */
public class ClassifierExample {
  /**
   * Entry point for the application. The arguments are ignored.
   */
  public static void main(String args[]) {
    JFrame f = new JFrame();
    f.addWindowListener(new WindowAdapter() {
      public void windowClosing(WindowEvent we) {
        System.exit(0);
      }
    });
    f.getContentPane().setLayout(new BorderLayout());
    final PointClassifier pc = new PointClassifier();
    f.getContentPane().add(BorderLayout.CENTER, pc);
    JPanel panel = new JPanel();
    panel.setLayout(new FlowLayout());
    ButtonGroup bGroup = new ButtonGroup();
    final JRadioButton rbPos = new JRadioButton("postive");
    bGroup.add(rbPos);
    final JRadioButton rbNeg = new JRadioButton("negative");
    bGroup.add(rbNeg);
    ActionListener addTypeAction = new ActionListener() {
      public void actionPerformed(ActionEvent ae) {
        //JRadioButton rb = (JRadioButton) ae.getSource();
        pc.setAddPos(rbPos.isSelected());
      }
    };
    rbPos.addActionListener(addTypeAction);
    panel.add(rbPos);
    rbNeg.addActionListener(addTypeAction);
    panel.add(rbNeg);
    ActionListener classifyAction = new ActionListener() {
      public void actionPerformed(ActionEvent ae) {
        pc.classify();
      }
    };
    JButton classifyB = new JButton("classify");
    classifyB.addActionListener(classifyAction);
    panel.add(classifyB);
    ActionListener clearAction = new ActionListener() {
      public void actionPerformed(ActionEvent ae) {
        pc.clear();
      }
    };
    JButton clearB = new JButton("clear");
    clearB.addActionListener(clearAction);
    panel.add(clearB);
    rbPos.setSelected(pc.getAddPos());
    rbNeg.setSelected(!pc.getAddPos());
    
    JComboBox kernelBox = new JComboBox();
    kernelBox.addItem("polynomeal");
    kernelBox.addItem("rbf");
    
    kernelBox.addItemListener(new ItemListener() {
      public void itemStateChanged(ItemEvent e) {
        if(e.getStateChange() == ItemEvent.SELECTED) {
          Object o = e.getItem();
          if(o.equals("polynomeal")) {
            pc.setKernel(PointClassifier.polyKernel);
          } else if(o.equals("rbf")) {
            pc.setKernel(PointClassifier.rbfKernel);
          }
        }
      }
    });
    panel.add(kernelBox);
    
    f.getContentPane().add(BorderLayout.NORTH, panel);
    f.setSize(400, 300);
    f.setVisible(true);
  }

  /**
   * An extention of JComponent that contains the points & encapsulates the
   * classifier.
   */
  public static class PointClassifier extends JComponent {
    // public kernels
    public static SVMKernel polyKernel;
    public static SVMKernel rbfKernel;
    public static SMOTrainer trainer;
    
    static {
      trainer = new SMOTrainer();
      trainer.setC(1.0E+7);
      trainer.setEpsilon(1.0E-9);
      
      SVMKernel k = new SVMKernel() {
        public double evaluate(Object a, Object b) {
          Point2D pa = (Point2D) a;
          Point2D pb = (Point2D) b;

          double dot = pa.getX() * pb.getX() + pa.getY() * pb.getY();
          return dot;
        }
      };

      PolynomialKernel pk = new PolynomialKernel();
      pk.setNestedKernel(k);
      pk.setOrder(2.0);
      pk.setConstant(1.0);
      pk.setMultiplier(0.0000001);
      
      RadialBaseKernel rb = new RadialBaseKernel();
      rb.setNestedKernel(k);
      rb.setWidth(10000.0);
      
      polyKernel = pk;
      rbfKernel = rb;
    }
    
    // private variables that should only be diddled by internal methods
    private SVMTarget target;
    private SVMClassifierModel model;

    {
      target = new SimpleSVMTarget();
      model = null;
    }

    // private variables containing state that may be diddled by beany methods
    private boolean addPos;
    private Shape posShape;
    private Shape negShape;
    private Paint svPaint;
    private Paint plainPaint;
    private Paint posPaint;
    private Paint negPaint;
    private SVMKernel kernel;

    /**
     * Set the kernel used for classification.
     *
     * @param kernel  the SVMKernel to use
     */
    public void setKernel(SVMKernel kernel) {
      firePropertyChange("kernel", this.kernel, kernel);
      this.kernel = kernel;
    }
    
    /**
     * Retrieve the currently used kernel
     *
     * @return the current value of the kernel.
     */
    public SVMKernel getKernel() {
      return this.kernel;
    }

    /**
     * Set a flag so that newly added points will be in the positive class or
     * negative class, depending on wether addPos is true or false respectively.
     *
     * @param addPos  boolean to flag which class to add new points to
     */
    public void setAddPos(boolean addPos) {
      firePropertyChange("addPos", this.addPos, addPos);
      this.addPos = addPos;
    }

    /**
     * Retrieve the current value of addPos.
     *
     * @return  true if new points will be added to the positive examples and
     *          false if they will be added to the negative examples.
     */
    public boolean getAddPos() {
      return addPos;
    }

    /**
     * Set the Shape to represent the positive points.
     * <p>
     * The shape should be positioned so that 0, 0 is the center or focus.
     *
     * @param posShape the Shape to use
     */
    public void setPosShape(Shape posShape) {
      firePropertyChange("posShape", this.posShape, posShape);
      this.posShape = posShape;
    }

    /**
     * Retrieve the shape used to represent positive points.
     *
     * @return the current positive Shape
     */
    public Shape getPosShape() {
      return posShape;
    }

    /**
     * Set the Shape to represent the negative points.
     * <p>
     * The shape should be positioned so that 0, 0 is the center or focus.
     *
     * @param negShape the Shape to use
     */
    public void setNegShape(Shape negShape) {
      firePropertyChange("negShape", this.negShape, negShape);
      this.negShape = negShape;
    }

    /**
     * Retrieve the shape used to represent negative points.
     *
     * @return the current negative Shape
     */
    public Shape getNegShape() {
      return negShape;
    }

    /**
     * Remove all points from the canvas, and discard any model.
     */
    public void clear() {
      target.clear();
      model = null;
      repaint();
    }

    /**
     * Learn a model from the current points.
     * <p>
     * This may take some time for complicated models.
     */
    public void classify() {
      new Thread() {
        public void run() {
          Cursor c = getCursor();
          setCursor(new Cursor(Cursor.WAIT_CURSOR));
          System.out.println("Training");
          model = trainer.trainModel(target, kernel, null);

          System.out.println("Threshold = " + model.getThreshold());
          for(Iterator i = model.items().iterator(); i.hasNext(); ) {
            Object item = i.next();
            System.out.println(item + "\t" +
                               target.getTarget(item) + "\t" +
                               model.getAlpha(item) + "\t" +
                               model.classify(item)
            );
          }

          PointClassifier.this.model = model;
          setCursor(c);
          repaint();
        }
      }.start();
    }

    /**
     * Make a new PointClassifier.
     * <p>
     * Hooks up the mouse listener & cursor.
     * Chooses default colors & Shapes.
     */
    public PointClassifier() {
      setCursor(new Cursor(Cursor.CROSSHAIR_CURSOR));
      addPos = true;
      setPosShape(new Rectangle2D.Double(-2.0, -2.0, 5.0, 5.0));
      setNegShape(new Ellipse2D.Double(-2.0, -2.0, 5.0, 5.0));
      setKernel(polyKernel);
      plainPaint = Color.black;
      svPaint = Color.green;
      posPaint = Color.red;
      negPaint = Color.blue;

      addMouseListener(new MouseAdapter() {
        public void mouseReleased(MouseEvent me) {
          Point p = me.getPoint();
          if(getAddPos()) {
            target.addItemTarget(p, +1.0);
          } else {
            target.addItemTarget(p, -1.0);
          }
          model = null;
          repaint();
        }
      });
    }
  
    /**
     * Renders this component to display the points, and if present, the
     * support vector machine.
     */
    public void paintComponent(Graphics g) {
      Graphics2D g2 = (Graphics2D) g;
      AffineTransform at = new AffineTransform();
      int i = 0;
      Rectangle r = g2.getClipBounds();
      int step = 3;

      if(model != null) {
        Rectangle rr = new Rectangle(r.x, r.y, step, step);
        Point p = new Point(r.x, r.y);
        for(int x = r.x; x < r.x + r.width; x+=step) {
          p.x = x;
          rr.x = x;
          for(int y = r.y; y < r.y + r.height; y+=step) {
            p.y = y;
            rr.y = y;
            double s = model.classify(p);
            if(s <= -1.0) {
              g2.setPaint(negPaint);
            } else if(s >= +1.0) {
              g2.setPaint(posPaint);
            } else {
              g2.setPaint(Color.white);
            }
            g2.fill(rr);
          }
        }
      }

      Set supportVectors = Collections.EMPTY_SET;
      if(model != null) {
        supportVectors = model.items();
      }
      for(Iterator it = target.items().iterator(); it.hasNext(); i++) {
        Point2D p = (Point2D) it.next();
        at.setToTranslation(p.getX(), p.getY());
        Shape glyph;
        if(target.getTarget(p) > 0) {
          glyph = getPosShape();
        } else {
          glyph = getNegShape();
        }
        Shape s = at.createTransformedShape(glyph);
        if(supportVectors.contains(p)) {
          g2.setPaint(svPaint);
        } else {
          g2.setPaint(plainPaint);
        }
        g2.draw(s);
      }
    }
  }
}
