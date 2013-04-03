package org.biojava3.structure.gui;

import java.awt.Color;
import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Rectangle;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.Box;
import javax.swing.JFrame;
import javax.swing.JMenuBar;
import javax.swing.JPanel;
import javax.swing.JTextField;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.gui.MenuCreator;

public class JmolViewerImpl implements StructureViewer {

    public static final String viewer = "org.jmol.api.JmolSimpleViewer";
    public static final String adapter = "org.jmol.api.JmolAdapter";
    public static final String smartAdapter = "org.jmol.adapter.smarter.SmarterJmolAdapter";
    Structure structure;
    JmolPanel jmolPanel;
    JFrame frame;

    public JmolViewerImpl() {

        frame = new JFrame();

        JMenuBar menu = MenuCreator.initMenu(frame, null,null);

        frame.setJMenuBar(menu);

        frame.addWindowListener(new WindowAdapter() {

            @Override
            public void windowClosing(WindowEvent e) {
                frame.dispose();
                //System.exit(0);
            }
        });

        Container contentPane = frame.getContentPane();

        Box vBox = Box.createVerticalBox();

        try {

            jmolPanel = new JmolPanel();

        } catch (ClassNotFoundException e) {
            e.printStackTrace();
            System.err.println("Could not find Jmol in classpath, please install first. http://www.jmol.org");
            return;
        }
        jmolPanel.setPreferredSize(new Dimension(500, 500));
        vBox.add(jmolPanel);


        JTextField field = new JTextField();

        field.setMaximumSize(new Dimension(Short.MAX_VALUE, 30));
        field.setText("enter RASMOL like command...");
//        RasmolCommandListener listener = new RasmolCommandListener(jmolPanel, field);

//        field.addActionListener(listener);
//        field.addMouseListener(listener);
//        field.addKeyListener(listener);
        vBox.add(field);

        contentPane.add(vBox);



        frame.pack();
        frame.setVisible(true);

    }

    public void setTitle(String label) {
        frame.setTitle(label);
        frame.repaint();
    }

    public void setStructure(Structure structure) {
        if (jmolPanel == null) {
            System.err.println("please install Jmol first");
            return;
        }

        setTitle(structure.getPDBCode());

        // actually this is very simple
        // just convert the structure to a PDB file

        String pdb = structure.toPDB();
        //System.out.println(s.isNmr());

        //System.out.println(pdb);
        // Jmol could also read the file directly from your file system
        //viewer.openFile("/Path/To/PDB/1tim.pdb");

        //System.out.println(pdb);
        jmolPanel.openStringInline(pdb);

        // send the PDB file to Jmol.
        // there are also other ways to interact with Jmol, e.g make it directly
        // access the biojava structure object, but they require more
        // code. See the SPICE code repository for how to do this.
    }

    public void clear() {
        // TODO Auto-generated method stub
    }

    public Color getColor() {
        // TODO Auto-generated method stub
        return null;
    }

    public Selection getSelection() {
        // TODO Auto-generated method stub
        return null;
    }

    public void repaint() {
        // TODO Auto-generated method stub
    }

    public void setColor(Color red) {
        // TODO Auto-generated method stub
    }

    public void setSelection(Selection selection) {
        // TODO Auto-generated method stub
    }

    public void setStyle(RenderStyle wireframe) {
        // TODO Auto-generated method stub
    }

    public void setZoom(int i) {
        // TODO Auto-generated method stub
    }
    @SuppressWarnings("rawtypes")
    static class JmolPanel extends JPanel {

        /**
         *
         */
        private static final long serialVersionUID = -3661941083797644242L;
        
		Class viewerC;
        
        Class adapterC;
        
        
		Class smartAdapterC;
        Object viewerO;
        Object adapterO;
        Method evalString;
        Method renderScreenImage;
        Method openStringInline;

        //JmolSimpleViewer viewer;
        //JmolAdapter adapter;
        @SuppressWarnings("unchecked")
		JmolPanel() throws ClassNotFoundException {

            try {
                viewerC = Class.forName(viewer);

                adapterC = Class.forName(adapter);
                smartAdapterC = Class.forName(smartAdapter);

                Method m = viewerC.getMethod("allocateSimpleViewer", new Class[]{Component.class, adapterC});

                Constructor constructor = smartAdapterC.getConstructor(new Class[]{});
                adapterO = constructor.newInstance(new Object[]{});

                //viewerC = JmolSimpleViewer.allocateSimpleViewer(this, adapter);
                viewerO = m.invoke(viewerC, this, adapterO);

                evalString = viewerC.getMethod("evalString", String.class);

                renderScreenImage = viewerC.getMethod("renderScreenImage",
                        new Class[]{Graphics.class, Dimension.class, Rectangle.class});

                openStringInline = viewerC.getMethod("openStringInline", new Class[]{String.class});

            } catch (InstantiationException ex) {
                Logger.getLogger(JmolViewerImpl.class.getName()).log(Level.SEVERE, null, ex);
            } catch (IllegalAccessException ex) {
                Logger.getLogger(JmolViewerImpl.class.getName()).log(Level.SEVERE, null, ex);
            } catch (IllegalArgumentException ex) {
                Logger.getLogger(JmolViewerImpl.class.getName()).log(Level.SEVERE, null, ex);
            } catch (InvocationTargetException ex) {
                Logger.getLogger(JmolViewerImpl.class.getName()).log(Level.SEVERE, null, ex);
            } catch (NoSuchMethodException e) {
                e.printStackTrace();
            }

            evalString("set scriptQueue on;");

        }

        
        public Class getViewer() {
            return viewerC;
        }

        public void evalString(String rasmolScript) {
            try {
                evalString.invoke(viewerO, rasmolScript);
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        public void openStringInline(String pdbFile) {
            try {
                openStringInline.invoke(viewerO, pdbFile);
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        public void executeCmd(String rasmolScript) {
            try {
                evalString.invoke(viewerO, rasmolScript);
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        final Dimension currentSize = new Dimension();
        final Rectangle rectClip = new Rectangle();

        public void paint(Graphics g) {
            getSize(currentSize);
            g.getClipBounds(rectClip);
            //viewer.renderScreenImage(g, currentSize, rectClip);

            try {
                renderScreenImage.invoke(viewerO, g, currentSize, rectClip);
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }
}
