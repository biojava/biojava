// $Id: MouseListener.java,v 1.4 2009/03/19 02:13:57 cmzmasek Exp $
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

import java.awt.Point;
import java.awt.event.InputEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionListener;

/*
 * @author Christian Zmasek
 */
final class MouseListener extends MouseAdapter implements MouseMotionListener {

    private final TreePanel _treepanel;
    private boolean         _being_dragged = false;
    private final Point     _click_point   = new Point();

    /**
     * Constructor.
     */
    MouseListener( final TreePanel tp ) {
        _treepanel = tp;
    }

    /**
     * Mouse clicked.
     */
    @Override
    public void mouseClicked( final MouseEvent e ) {
        _click_point.setLocation( e.getX(), e.getY() );
        _treepanel.mouseClicked( e );
    }

    @Override
    public void mouseDragged( final MouseEvent e ) {
        if ( ( e.getModifiersEx() == InputEvent.BUTTON1_DOWN_MASK )
                || ( e.getModifiersEx() == InputEvent.BUTTON3_DOWN_MASK ) ) {
            if ( !_treepanel.inOvRectangle( e ) ) {
                if ( !_being_dragged ) {
                    _being_dragged = true;
                    _treepanel.setLastMouseDragPointX( e.getX() );
                    _treepanel.setLastMouseDragPointY( e.getY() );
                }
                _treepanel.mouseDragInBrowserPanel( e );
            }
            else {
                if ( !_being_dragged ) {
                    _being_dragged = true;
                    _treepanel.setLastMouseDragPointX( e.getX() );
                    _treepanel.setLastMouseDragPointY( e.getY() );
                }
                _treepanel.mouseDragInOvRectangle( e );
            }
        }
    }

    @Override
    public void mouseMoved( final MouseEvent e ) {
        _treepanel.mouseMoved( e );
    }

    @Override
    public void mousePressed( final MouseEvent e ) {
        //TODO is this a good idea? It is certainly not NEEDED.
        if ( e.getModifiersEx() == InputEvent.BUTTON1_DOWN_MASK ) {
            if ( !_being_dragged ) {
                _being_dragged = true;
                _treepanel.setLastMouseDragPointX( e.getX() );
                _treepanel.setLastMouseDragPointY( e.getY() );
            }
            if ( !_treepanel.inOvRectangle( e ) ) {
                _treepanel.mouseDragInBrowserPanel( e );
            }
            else {
                _treepanel.mouseDragInOvRectangle( e );
            }
        }
    }

    @Override
    public void mouseReleased( final MouseEvent e ) {
        if ( _being_dragged ) {
            _being_dragged = false;
        }
        _treepanel.mouseReleasedInBrowserPanel( e );
    }
}