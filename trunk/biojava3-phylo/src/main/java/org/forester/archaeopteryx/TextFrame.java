// $Id: TextFrame.java,v 1.4 2009/10/26 23:29:39 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// Copyright (C) 2000-2001 Washington University School of Medicine
// and Howard Hughes Medical Institute
// Copyright (C) 2003-2007 Ethalinda K.S. Cannon
// All rights reserved
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: cmzmasek@yahoo.com
// WWW: www.phylosoft.org/forester

package org.forester.archaeopteryx;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.ClipboardOwner;
import java.awt.datatransfer.StringSelection;
import java.awt.datatransfer.Transferable;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

final class TextFrame extends JFrame implements ActionListener, ClipboardOwner {

    /**
     * 
     */
    private static final long serialVersionUID = -5012834229705518363L;
    private static Color      ta_text_color    = new Color( 0, 0, 0 ),
            ta_background_color = new Color( 240, 240, 240 ), background_color = new Color( 215, 215, 215 ),
            button_background_color = new Color( 215, 215, 215 ), button_text_color = new Color( 0, 0, 0 );
    private final static Font button_font      = new Font( "Helvetica", Font.PLAIN, 10 ),
            ta_font = new Font( "Helvetica", Font.PLAIN, 10 );
    private boolean           can_use_clipboard;
    private final String      text;
    private final JTextArea   jtextarea;
    private final JButton     close_button;
    private JButton           copy_button;
    private final JPanel      buttonjpanel;
    private final Container   contentpane;

    private TextFrame( final String s ) {
        // first things first
        setTitle( Constants.PRG_NAME );
        text = s;
        // check to see if we have permission to use the clipboard:
        can_use_clipboard = true;
        final SecurityManager sm = System.getSecurityManager();
        if ( sm != null ) {
            try {
                sm.checkSystemClipboardAccess();
            }
            catch ( final Exception e ) {
                //nope!
                can_use_clipboard = false;
            }
        }
        // set up the frame
        setBackground( background_color );
        buttonjpanel = new JPanel();
        buttonjpanel.setBackground( background_color );
        close_button = new JButton( "          Close          " );
        close_button.setBackground( button_background_color );
        close_button.setForeground( button_text_color );
        close_button.setFont( button_font );
        close_button.addActionListener( this );
        buttonjpanel.add( close_button );
        if ( can_use_clipboard ) {
            copy_button = new JButton( "Copy to clipboard" );
            copy_button.setBackground( button_background_color );
            copy_button.setForeground( button_text_color );
            copy_button.setFont( button_font );
            copy_button.addActionListener( this );
            buttonjpanel.add( copy_button );
        }
        contentpane = getContentPane();
        contentpane.setLayout( new BorderLayout() );
        jtextarea = new JTextArea( text );
        jtextarea.setBackground( ta_background_color );
        jtextarea.setForeground( ta_text_color );
        jtextarea.setFont( ta_font );
        jtextarea.setEditable( false );
        jtextarea.setWrapStyleWord( true );
        jtextarea.setLineWrap( true );
        contentpane.add( new JScrollPane( jtextarea ), BorderLayout.CENTER );
        buttonjpanel.setLayout( new FlowLayout( FlowLayout.CENTER, 20, 5 ) );
        contentpane.add( buttonjpanel, BorderLayout.SOUTH );
        setSize( 500, 400 );
        addWindowListener( new WindowAdapter() {

            @Override
            public void windowClosing( final WindowEvent e ) {
                close();
            }
        } );
        setVisible( true );
    }

    public void actionPerformed( final ActionEvent e ) {
        final Object o = e.getSource();
        if ( o == close_button ) {
            close();
        }
        else if ( o == copy_button ) {
            copy();
        }
    }

    void close() {
        setVisible( false );
        dispose();
    }

    private void copy() {
        if ( !can_use_clipboard ) {
            // can't do this!
            return;
        }
        final Clipboard sys_clipboard = getToolkit().getSystemClipboard();
        final StringSelection contents = new StringSelection( jtextarea.getText() );
        sys_clipboard.setContents( contents, this );
    }

    public void lostOwnership( final Clipboard clipboard, final Transferable contents ) {
    }

    static TextFrame instantiate( final String s ) {
        return new TextFrame( s );
    }
}
