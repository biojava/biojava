// $Id: NodeFrame.java,v 1.7 2009/10/27 00:26:46 cmzmasek Exp $
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
import java.awt.Container;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.util.ForesterUtil;

final class NodeFrame extends javax.swing.JFrame {

    private static final long serialVersionUID = -6943510233968557246L;
    private final TreePanel   _reepanel;
    private int               _index           = -1;

    NodeFrame( final PhylogenyNode n, final Phylogeny tree, final TreePanel tp, final int x ) {
        super( "Node " + ( ForesterUtil.isEmpty( n.getNodeName() ) ? n.getNodeId() : n.getNodeName() ) );
        _reepanel = tp;
        setSize( Constants.NODE_FRAME_SIZE );
        _index = x;
        final Container contentPane = getContentPane();
        final NodePanel nodepanel = new NodePanel( n );
        contentPane.add( nodepanel, BorderLayout.CENTER );
        addWindowListener( new WindowAdapter() {

            @Override
            public void windowClosing( final WindowEvent e ) {
                remove(); // to release slot in array
                dispose();
            }
        } );
        setResizable( false );
        nodepanel.setVisible( true );
        setVisible( true );
    }

    NodeFrame( final PhylogenyNode n, final Phylogeny tree, final TreePanel tp, final int x, final String dummy ) {
        super( "Editable Node " + ( ForesterUtil.isEmpty( n.getNodeName() ) ? n.getNodeId() : n.getNodeName() ) );
        _reepanel = tp;
        setSize( Constants.NODE_FRAME_SIZE );
        _index = x;
        final Container contentPane = getContentPane();
        final NodeEditPanel nodepanel = new NodeEditPanel( n, tp );
        contentPane.add( nodepanel, BorderLayout.CENTER );
        addWindowListener( new WindowAdapter() {

            @Override
            public void windowClosing( final WindowEvent e ) {
                try {
                    nodepanel.writeAll();
                }
                catch ( final Exception ex ) {
                    // Do nothing.
                }
                remove(); // to release slot in array
                dispose();
            }
        } );
        setResizable( false );
        nodepanel.setVisible( true );
        setVisible( true );
    }

    TreePanel getTreePanel() {
        return _reepanel;
    }

    void remove() {
        if ( _index > -1 ) {
            _reepanel.removeEditNodeFrame( _index ); // to release slot in array
        }
    }
}
