/**
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
 * Created on Sep 14, 2011
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.bio.structure.align.gui.autosuggest;

import java.awt.Font;
import java.awt.Frame;
import java.awt.IllegalComponentStateException;
import java.awt.Point;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.util.Iterator;
import java.util.Vector;
import java.util.concurrent.atomic.AtomicBoolean;

import javax.swing.JDialog;
import javax.swing.JList;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;
import javax.swing.SwingWorker;




/** A JTextField that can make suggestions for auto-complete.
 * 
 * @author Andreas Prlic
 *
 */
public class JAutoSuggest extends JTextField{

	/**
	 * 
	 */
	private static final long serialVersionUID = 8591734727984365156L;

	private static final String DEFAULT_TEXT= "Please enter text ...";
	
	String defaultText;
	private JDialog dialog;	
	private Point location;
	private JList list;
	
	private Vector<String> suggestions;
	
	/** last word that was entered by user */
	private String lastWord ;
	
	AutoSuggestProvider autoSuggestProvider;

	Font regular;
	Font busy;
	
	SuggestionFetcher matcher;
	
	public JAutoSuggest(){
		super();
		
		init();
	}
	
	public JAutoSuggest(int size){
		super(size);
		init();
	}
	
	public JAutoSuggest(Frame owner){
		owner.addComponentListener(new ComponentListener() {
			@Override
			public void componentShown(ComponentEvent e) {
				updateLocation();
			}

			@Override
			public void componentResized(ComponentEvent e) {
				updateLocation();
			}

			@Override
			public void componentMoved(ComponentEvent e) {
				updateLocation();
			}

			@Override
			public void componentHidden(ComponentEvent e) {
				updateLocation();
			}
		});
		owner.addWindowListener(new WindowListener() {
			@Override
			public void windowOpened(WindowEvent e) {
			}

			@Override
			public void windowIconified(WindowEvent e) {
				dialog.setVisible(false);
			}

			@Override
			public void windowDeiconified(WindowEvent e) {
			}

			@Override
			public void windowDeactivated(WindowEvent e) {
			}

			@Override
			public void windowClosing(WindowEvent e) {
				dialog.dispose();
			}

			@Override
			public void windowClosed(WindowEvent e) {
				dialog.dispose();
			}

			@Override
			public void windowActivated(WindowEvent e) {
			}
		});
		addFocusListener(new FocusListener() {
			@Override
			public void focusLost(FocusEvent e) {
				dialog.setVisible(false);
				
				if (getText().equals("") && e.getOppositeComponent() != null && e.getOppositeComponent().getName() != null) {
					if (!e.getOppositeComponent().getName().equals("suggestFieldDropdownButton")) {
						setText(defaultText);
					}
				} else if (getText().equals("")) {
					setText(defaultText);
				}
			}

			@Override
			public void focusGained(FocusEvent e) {
				if (getText().equals(defaultText)) {
					setText("");
				}
				
				showSuggest();
			}
		});
		
		
		
		
		init();
		
		// set dialog owner...
		//dialog.setD
		
	}
	
	
	private void initSuggestionList() {
		list = new JList();
		list.addMouseListener(new MouseListener() {
			private int selected;

			@Override
			public void mousePressed(MouseEvent e) {
			}

			@Override
			public void mouseReleased(MouseEvent e) {
				if (selected == list.getSelectedIndex()) {
					// provide double-click for selecting a suggestion
					setText((String) list.getSelectedValue());					

					dialog.setVisible(false);
				}
				selected = list.getSelectedIndex();
			}

			@Override
			public void mouseExited(MouseEvent e) {
			}

			@Override
			public void mouseEntered(MouseEvent e) {
			}

			@Override
			public void mouseClicked(MouseEvent e) {
			}
		});
		dialog.add(new JScrollPane(list, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
				JScrollPane.HORIZONTAL_SCROLLBAR_NEVER));
		dialog.pack();
		addKeyListener(new KeyListener() {
			@Override
			public void keyTyped(KeyEvent e) {
			}

			@Override
			public void keyPressed(KeyEvent e) {
				updateLocation();
			}

			@Override
			public void keyReleased(KeyEvent e) {
				if (e.getKeyCode() == KeyEvent.VK_ESCAPE) {
					dialog.setVisible(false);
					return;
				} else if (e.getKeyCode() == KeyEvent.VK_DOWN) {
					if (dialog.isVisible()) {
						list.setSelectedIndex(list.getSelectedIndex() + 1);
						list.ensureIndexIsVisible(list.getSelectedIndex() + 1);
						return;
					} else {
						showSuggest();
					}
				} else if (e.getKeyCode() == KeyEvent.VK_UP) {
					list.setSelectedIndex(list.getSelectedIndex() - 1);
					list.ensureIndexIsVisible(list.getSelectedIndex() - 1);
					return;
				} else if (e.getKeyCode() == KeyEvent.VK_ENTER
						& list.getSelectedIndex() != -1 & suggestions.size() > 0) {
					setText((String) list.getSelectedValue());
					
					
					dialog.setVisible(false);
					return;
				}
				showSuggest();
			}
		});
		
	}

	private void init(){
		autoSuggestProvider =  new DefaultAutoSuggestProvider();
		lastWord = "";
		regular = getFont();
		busy = new Font(getFont().getName(), Font.ITALIC, getFont().getSize());
		suggestions = new Vector<String>();
		defaultText = DEFAULT_TEXT;
		
		
		dialog = new JDialog();
		dialog.setUndecorated(true);
		dialog.setFocusableWindowState(false);
		dialog.setFocusable(false);
		
		initSuggestionList();
	}
	
	public String getDefaultText() {
		return defaultText;
	}

	public void setDefaultText(String defaultText) {
		this.defaultText = defaultText;
	}

	public AutoSuggestProvider getAutoSuggestProvider() {
		return autoSuggestProvider;
	}

	public void setAutoSuggestProvider(AutoSuggestProvider autoSuggestProvider) {
		this.autoSuggestProvider = autoSuggestProvider;
	}

	/**
	 * Force the suggestions to be displayed (Useful for buttons
	 * e.g. for using JSuggestionField like a ComboBox)
	 */
	public void showSuggest() {
		
		lastWord = getText();
		if ( lastWord != null) {
		
		}
			//autoSuggestProvider.getSuggestion(lastWord);
		
		
		
		if (!getText().toLowerCase().contains(lastWord.toLowerCase())) {
			suggestions.clear();
		}
		
		if (matcher != null) {
			matcher.setStop();
		}
		matcher = new SuggestionFetcher();

		//SwingUtilities.invokeLater(matcher);
		matcher.execute();
		lastWord = getText();
		updateLocation();
	}

	/**
	 * Force the suggestions to be hidden (Useful for buttons, e.g. to use
	 * JSuggestionField like a ComboBox)
	 */
	public void hideSuggest() {
		dialog.setVisible(false);
	}

	/**
	 * @return boolean Visibility of the suggestion window
	 */
	public boolean isSuggestVisible() {
		return dialog.isVisible();
	}

	
	
	/**
	 * Place the suggestion window under the JTextField.
	 */
	private void updateLocation() {
		try {
			location = getLocationOnScreen();
			location.y += getHeight();
			dialog.setLocation(location);
		} catch (IllegalComponentStateException e) {
			return; // might happen on window creation
		}
	}
	
	
	/** fetch suggestions from SuggestionProvider
	 * 
	 * 
	 *
	 */
	private class SuggestionFetcher extends SwingWorker<String, Object> {
		/** flag used to stop the thread */
		private AtomicBoolean stop = new AtomicBoolean(false);
		
		String previousWord;
		/**
		 * Standard run method used in threads
		 * responsible for the actual search
		 */
		@Override
		public String doInBackground() {
			try {
				setFont(busy);
				String userInput = getText();
				if ( userInput == null || userInput.equals(""))
					return "";
				
				if ( previousWord != null){
					if ( userInput.equals(previousWord))
						return "";
				}
				previousWord = userInput;
			
				suggestions = autoSuggestProvider.getSuggestion(userInput);
							
				setFont(regular);
				
				if (suggestions.size() > 0) {
					list.setListData(suggestions);
					list.setSelectedIndex(0);
					list.ensureIndexIsVisible(0);
					dialog.setVisible(true);
				} else {
					dialog.setVisible(false);
				}
			} catch (Exception e) {
				//e.printStackTrace();
				// ignore...
				
			}
			return "Done.";
		}

		public void setStop(){
			stop.set(true);
		}
		
		
	}
}
