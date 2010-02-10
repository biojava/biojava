// $Id: ArchaeopteryxA.java,v 1.6 2009/11/20 22:22:10 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
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

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.KeyboardFocusManager;
import java.io.File;
import java.net.URL;

import javax.swing.JApplet;
import javax.swing.UIManager;

import org.forester.phylogeny.Phylogeny;
import org.forester.util.ForesterUtil;

public class ArchaeopteryxA extends JApplet {

    private static final long  serialVersionUID    = 2314899014580484146L;
    private final static Color background_color    = new Color( 0, 0, 0 );
    private final static Color font_color          = new Color( 0, 255, 0 );
    private final static Color ex_background_color = new Color( 0, 0, 0 );
    private final static Color ex_font_color       = new Color( 255, 0, 0 );
    private final static Font  font                = new Font( Configuration.getDefaultFontFamilyName(), Font.BOLD, 9 );
    private MainFrameApplet    _mainframe_applet;
    private String             _url_string         = "";
    private String             _message_1          = "";
    private String             _message_2          = "";
    public final static String NAME                = "ArchaeopteryxA";

    @Override
    public void destroy() {
        Util.printAppletMessage( NAME, "going to be destroyed" );
        if ( getMainFrameApplet() != null ) {
            getMainFrameApplet().close();
        }
    }

    private MainFrameApplet getMainFrameApplet() {
        return _mainframe_applet;
    }

    private String getMessage1() {
        return _message_1;
    }

    private String getMessage2() {
        return _message_2;
    }

    public String getUrlString() {
        return _url_string;
    }

    @Override
    public void init() {
        boolean has_exception = false;
        setName( NAME );
        setUrlString( getParameter( Constants.APPLET_PARAM_NAME_FOR_URL_OF_TREE_TO_LOAD ) );
        Util.printAppletMessage( NAME, "URL of phylogenies to load: \"" + getUrlString() + "\"" );
        setBackground( background_color );
        setForeground( font_color );
        setFont( font );
        repaint();
        String s = null;
        try {
            s = System.getProperty( "java.version" );
        }
        catch ( final Exception e ) {
            ForesterUtil.printWarningMessage( NAME, "minor error: " + e.getLocalizedMessage() );
        }
        if ( ( s != null ) && ( s.length() > 0 ) ) {
            setMessage2( "[Your Java version: " + s + "]" );
            repaint();
        }
        final String config_filename = getParameter( Constants.APPLET_PARAM_NAME_FOR_CONFIG_FILE_URL );
        Util.printAppletMessage( NAME, "URL for configuration file is: " + config_filename );
        final Configuration configuration = new Configuration( config_filename, true, true );
        try {
            if ( configuration.isUseNativeUI() ) {
                UIManager.setLookAndFeel( UIManager.getSystemLookAndFeelClassName() );
            }
            else {
                UIManager.setLookAndFeel( UIManager.getCrossPlatformLookAndFeelClassName() );
            }
            setVisible( false );
            _mainframe_applet = new MainFrameApplet( this, configuration );
            URL url = null;
            url = new URL( getUrlString() );
            final Phylogeny[] phys = Util.readPhylogeniesFromUrl( url, configuration.isValidatePhyloXmlAgainstSchema() );
            Util.addPhylogeniesToTabs( phys, new File( url.getFile() ).getName(), getUrlString(), getMainFrameApplet()
                    .getConfiguration(), getMainFrameApplet().getMainPanel() );
            getMainFrameApplet().getMainPanel().getControlPanel().showWholeAll();
            getMainFrameApplet().getMainPanel().getControlPanel().showWhole();
            setVisible( true );
        }
        catch ( final Exception e ) {
            ForesterUtil.printErrorMessage( NAME, e.toString() );
            setBackground( ex_background_color );
            setForeground( ex_font_color );
            has_exception = true;
            setMessage1( "Exception: " + e );
            e.printStackTrace();
            repaint();
        }
        if ( !has_exception ) {
            setMessage1( NAME + " is now ready!" );
            repaint();
            Util.printAppletMessage( NAME, "successfully initialized" );
        }
        KeyboardFocusManager.getCurrentKeyboardFocusManager().clearGlobalFocusOwner();
        getMainFrameApplet().requestFocus();
        getMainFrameApplet().requestFocusInWindow();
        getMainFrameApplet().requestFocus();
    }

    /**
     * Prints message when initialization is finished. Called automatically.
     * 
     * @param g
     *            Graphics
     */
    @Override
    public void paint( final Graphics g ) {
        g.setColor( background_color );
        g.fillRect( 0, 0, 300, 60 );
        g.setColor( font_color );
        g.drawString( getMessage2(), 10, 20 );
        g.drawString( getMessage1(), 10, 40 );
    }

    private void setMessage1( final String message_1 ) {
        _message_1 = message_1;
    }

    private void setMessage2( final String message_2 ) {
        _message_2 = message_2;
    }

    private void setUrlString( final String url_string ) {
        _url_string = url_string;
    }

    @Override
    public void start() {
        getMainFrameApplet().getMainPanel().validate();
        getMainFrameApplet().requestFocus();
        getMainFrameApplet().requestFocusInWindow();
        getMainFrameApplet().requestFocus();
        Util.printAppletMessage( NAME, "started" );
    }
}
