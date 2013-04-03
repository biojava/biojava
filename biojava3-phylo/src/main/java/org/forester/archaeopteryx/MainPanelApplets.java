// $Id: MainPanelApplets.java,v 1.6 2009/03/08 05:56:37 cmzmasek Exp $
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
import java.util.ArrayList;

import javax.swing.JApplet;

final class MainPanelApplets extends MainPanel {

    private static final long serialVersionUID = -7142615479464963140L;
    private final JApplet     _applet;

    public MainPanelApplets( final Configuration configuration, final ArchaeopteryxE em_applet ) {
        if ( configuration == null ) {
            throw new IllegalArgumentException( "configuration is null" );
        }
        addComponentListener( this );
        _configuration = configuration;
        _mainframe = null;
        _treepanels = new ArrayList<TreePanel>();
        _applet = em_applet;
        initialize();
        _control_panel = new ControlPanel( this, configuration );
        if ( !configuration.isHideControlPanelAndMenubar() ) {
            add( _control_panel, BorderLayout.WEST );
        }
        setupTreeGraphic( configuration, getControlPanel() );
    }

    public MainPanelApplets( final Configuration configuration, final MainFrameApplet aaf ) {
        if ( configuration == null ) {
            throw new IllegalArgumentException( "configuration is null" );
        }
        addComponentListener( this );
        _configuration = configuration;
        _mainframe = aaf;
        _treepanels = new ArrayList<TreePanel>();
        _applet = aaf.getApplet();
        initialize();
        _control_panel = new ControlPanel( this, configuration );
        add( _control_panel, BorderLayout.WEST );
        setupTreeGraphic( configuration, getControlPanel() );
    }

    JApplet getApplet() {
        return _applet;
    }

    MainFrameApplet getAppletFrame() {
        return ( MainFrameApplet ) _mainframe;
    }

    @Override
    Options getOptions() {
        if ( _mainframe != null ) {
            return _mainframe.getOptions();
        }
        else {
            return ( ( ArchaeopteryxE ) _applet ).getOptions();
        }
    }
}