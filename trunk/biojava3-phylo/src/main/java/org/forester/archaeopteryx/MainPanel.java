// $Id: MainPanel.java,v 1.19 2009/11/17 03:46:45 cmzmasek Exp $
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
import java.awt.Dimension;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.SwingConstants;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.forester.archaeopteryx.phylogeny.data.RenderableDomainArchitecture;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.util.ForesterUtil;

public class MainPanel extends JPanel implements ComponentListener {

    private static final long  serialVersionUID = -2682765312661416435L;
    MainFrame                  _mainframe;
    List<TreePanel>            _treepanels;
    ControlPanel               _control_panel;
    private List<JScrollPane>  _treegraphic_scroll_panes;
    private List<JPanel>       _treegraphic_scroll_pane_panels;
    Configuration              _configuration;
    private JTabbedPane        _tabbed_pane;
    private TreeColorSet       _colorset;
    private TreeFontSet        _fontset;
    private Phylogeny          _cut_or_copied_tree;
    private Set<PhylogenyNode> _copied_and_pasted_nodes;

    MainPanel() {
    }

    MainPanel( final Configuration configuration, final MainFrame parent ) {
        if ( configuration == null ) {
            throw new IllegalArgumentException( "configuration is null" );
        }
        addComponentListener( this );
        _configuration = configuration;
        _mainframe = parent;
        _treepanels = new ArrayList<TreePanel>();
        initialize();
        _control_panel = new ControlPanel( this, configuration );
        add( _control_panel, BorderLayout.WEST );
        setupTreeGraphic( configuration, getControlPanel() );
        getControlPanel().showWhole();
    }

    void addPhylogenyInNewTab( final Phylogeny phy,
                               final Configuration config,
                               final String default_name,
                               final String full_path ) {
        final TreePanel treepanel = new TreePanel( phy, config, this );
        getControlPanel().phylogenyAdded( config );
        treepanel.setControlPanel( getControlPanel() );
        _treepanels.add( treepanel );
        String name = "";
        if ( !ForesterUtil.isEmpty( phy.getName() ) ) {
            name = phy.getName();
        }
        else if ( phy.getIdentifier() != null ) {
            name = phy.getIdentifier().toString();
        }
        else if ( !ForesterUtil.isEmpty( default_name ) ) {
            name = default_name;
        }
        else {
            name = "[" + ( getTabbedPane().getTabCount() + 1 ) + "]";
        }
        final JScrollPane treegraphic_scroll_pane = new JScrollPane( treepanel );
        treegraphic_scroll_pane.getHorizontalScrollBar().addAdjustmentListener( new AdjustmentListener() {

            @Override
            public void adjustmentValueChanged( final AdjustmentEvent e ) {
                if ( treepanel.isOvOn() || getOptions().isShowScale() ) {
                    treepanel.repaint();
                }
            }
        } );
        treegraphic_scroll_pane.getVerticalScrollBar().addAdjustmentListener( new AdjustmentListener() {

            @Override
            public void adjustmentValueChanged( final AdjustmentEvent e ) {
                if ( treepanel.isOvOn() || getOptions().isShowScale() ) {
                    treepanel.repaint();
                }
            }
        } );
        treegraphic_scroll_pane.getHorizontalScrollBar().setUnitIncrement( 10 );
        treegraphic_scroll_pane.getHorizontalScrollBar().setBlockIncrement( 200 );
        treegraphic_scroll_pane.getVerticalScrollBar().setUnitIncrement( 10 );
        treegraphic_scroll_pane.getVerticalScrollBar().setBlockIncrement( 200 );
        final JPanel treegraphic_scroll_pane_panel = new JPanel();
        treegraphic_scroll_pane_panel.setLayout( new BorderLayout() );
        treegraphic_scroll_pane_panel.add( treegraphic_scroll_pane, BorderLayout.CENTER );
        _treegraphic_scroll_pane_panels.add( treegraphic_scroll_pane_panel );
        _treegraphic_scroll_panes.add( treegraphic_scroll_pane );
        getTabbedPane().addTab( name,
                                null,
                                treegraphic_scroll_pane_panel,
                                Util.createDescriptionForTab( phy, full_path ) );
        getTabbedPane().setSelectedIndex( getTabbedPane().getTabCount() - 1 );
        getControlPanel().showWhole();
    }

    void addPhylogenyInPanel( final Phylogeny phy, final Configuration config ) {
        final TreePanel treepanel = new TreePanel( phy, config, this );
        getControlPanel().phylogenyAdded( config );
        treepanel.setControlPanel( getControlPanel() );
        _treepanels.add( treepanel );
        final JScrollPane treegraphic_scroll_pane = new JScrollPane( treepanel );
        treegraphic_scroll_pane.getHorizontalScrollBar().setUnitIncrement( 20 );
        treegraphic_scroll_pane.getHorizontalScrollBar().setBlockIncrement( 50 );
        treegraphic_scroll_pane.getVerticalScrollBar().setUnitIncrement( 20 );
        treegraphic_scroll_pane.getVerticalScrollBar().setBlockIncrement( 50 );
        final JPanel treegraphic_scroll_pane_panel = new JPanel();
        treegraphic_scroll_pane_panel.setLayout( new BorderLayout() );
        treegraphic_scroll_pane_panel.add( treegraphic_scroll_pane, BorderLayout.CENTER );
        _treegraphic_scroll_pane_panels.add( treegraphic_scroll_pane_panel );
        _treegraphic_scroll_panes.add( treegraphic_scroll_pane );
        add( treegraphic_scroll_pane_panel, BorderLayout.CENTER );
    }

    void adjustJScrollPane() {
        if ( getTabbedPane() != null ) {
            getCurrentScrollPanePanel().remove( getCurrentScrollPane() );
            getCurrentScrollPanePanel().add( getCurrentScrollPane(), BorderLayout.CENTER );
        }
        getCurrentScrollPane().revalidate();
    }

    void closeCurrentPane() {
        final int index = getCurrentTabIndex();
        if ( ( index >= 0 ) && ( getTabbedPane().getTabCount() > 0 ) ) {
            getTabbedPane().remove( index );
            getTreePanels().remove( index );
            _treegraphic_scroll_panes.remove( index );
            _treegraphic_scroll_pane_panels.remove( index );
            getControlPanel().phylogenyRemoved( index );
        }
    }

    @Override
    public void componentHidden( final ComponentEvent e ) {
        // Do nothing.
    }

    @Override
    public void componentMoved( final ComponentEvent e ) {
        // Do nothing.
    }

    @Override
    public void componentResized( final ComponentEvent e ) {
        if ( getCurrentTreePanel() != null ) {
            getCurrentTreePanel().updateOvSettings();
            getCurrentTreePanel().updateOvSizes();
        }
    }

    @Override
    public void componentShown( final ComponentEvent e ) {
        // Do nothing.
    }

    private Configuration getConfiguration() {
        return _configuration;
    }

    ControlPanel getControlPanel() {
        return _control_panel;
    }

    public Set<PhylogenyNode> getCopiedAndPastedNodes() {
        return _copied_and_pasted_nodes;
    }

    Phylogeny getCurrentPhylogeny() {
        if ( getCurrentTreePanel() == null ) {
            return null;
        }
        return getCurrentTreePanel().getPhylogeny();
    }

    JScrollPane getCurrentScrollPane() {
        if ( _treegraphic_scroll_panes.size() > 0 ) {
            final int selected = _tabbed_pane.getSelectedIndex();
            if ( selected >= 0 ) {
                return _treegraphic_scroll_panes.get( selected );
            }
            else {
                return _treegraphic_scroll_panes.get( 0 );
            }
        }
        else {
            return null;
        }
    }

    JPanel getCurrentScrollPanePanel() {
        final int selected = _tabbed_pane.getSelectedIndex();
        if ( selected >= 0 ) {
            return _treegraphic_scroll_pane_panels.get( selected );
        }
        else {
            return _treegraphic_scroll_pane_panels.get( 0 );
        }
    }

    int getCurrentTabIndex() {
        final int selected = _tabbed_pane.getSelectedIndex();
        if ( selected >= 0 ) {
            return selected;
        }
        else {
            return 0;
        }
    }

    TreePanel getCurrentTreePanel() {
        final int selected = getTabbedPane().getSelectedIndex();
        if ( selected >= 0 ) {
            return _treepanels.get( selected );
        }
        else {
            if ( _treepanels.size() == 1 ) {
                return _treepanels.get( 0 );
            }
            else {
                return null;
            }
        }
    }

    Phylogeny getCutOrCopiedTree() {
        return _cut_or_copied_tree;
    }

    MainFrame getMainFrame() {
        return _mainframe;
    }

    Options getOptions() {
        return _mainframe.getOptions();
    }

    Phylogeny getPhylogeny( final int index ) {
        if ( getCurrentTreePanel() == null ) {
            return null;
        }
        return _treepanels.get( index ).getPhylogeny();
    }

    Dimension getSizeOfViewport() {
        return getCurrentScrollPane().getViewport().getExtentSize();
    }

    JTabbedPane getTabbedPane() {
        return _tabbed_pane;
    }

    TreeColorSet getTreeColorSet() {
        return _colorset;
    }

    public TreeFontSet getTreeFontSet() {
        return _fontset;
    }

    List<TreePanel> getTreePanels() {
        return _treepanels;
    }

    void initialize() {
        if ( !getConfiguration().isUseNativeUI() ) {
            setBackground( Constants.GUI_BACKGROUND_DEFAULT );
        }
        setTreeFontSet( new TreeFontSet( this ) );
        getTreeFontSet().setBaseFont( getOptions().getBaseFont() );
        setLayout( new BorderLayout() );
        setTreeColorSet( TreeColorSet.createInstance( getConfiguration() ) );
        _treegraphic_scroll_panes = new ArrayList<JScrollPane>();
        _treegraphic_scroll_pane_panels = new ArrayList<JPanel>();
        _tabbed_pane = new JTabbedPane( SwingConstants.TOP );
        if ( !getConfiguration().isUseNativeUI() ) {
            _tabbed_pane.setBackground( Constants.GUI_BACKGROUND_DEFAULT );
            _tabbed_pane.setForeground( Constants.GUI_BACKGROUND_DEFAULT );
        }
        _tabbed_pane.addChangeListener( new ChangeListener() {

            // This method is called whenever the selected tab changes
            public void stateChanged( final ChangeEvent evt ) {
                final JTabbedPane pane = ( JTabbedPane ) evt.getSource();
                getControlPanel().tabChanged();
                // Get current tab
                final int sel = pane.getSelectedIndex();
                if ( sel >= 0 ) {
                    if ( !getConfiguration().isUseNativeUI() ) {
                        if ( _tabbed_pane.getTabCount() > 0 ) {
                            _tabbed_pane.setForegroundAt( sel, Constants.TAB_LABEL_FOREGROUND_COLOR_SELECTED );
                            for( int i = 0; i < _tabbed_pane.getTabCount(); ++i ) {
                                if ( i != sel ) {
                                    _tabbed_pane.setBackgroundAt( i, Constants.TAB_LABEL_BACKGROUND_COLOR_NOT_SELECTED );
                                    _tabbed_pane.setForegroundAt( i, Constants.TAB_LABEL_FOREGROUND_COLOR_NOT_SELECTED );
                                }
                            }
                        }
                    }
                }
            }
        } );
        if ( !getConfiguration().isUseNativeUI() ) {
            _tabbed_pane.setFont( ControlPanel.jcb_font );
        }
        _tabbed_pane.setTabLayoutPolicy( JTabbedPane.SCROLL_TAB_LAYOUT );
        add( _tabbed_pane, BorderLayout.CENTER );
    }

    public void setArrowCursor() {
        setCursor( TreePanel.ARROW_CURSOR );
        repaint();
    }

    public void setCopiedAndPastedNodes( final Set<PhylogenyNode> copied_and_pasted_nodes ) {
        _copied_and_pasted_nodes = copied_and_pasted_nodes;
    }

    void setCutOrCopiedTree( final Phylogeny cut_or_copied_tree ) {
        _cut_or_copied_tree = cut_or_copied_tree;
    }

    void setTreeColorSet( final TreeColorSet colorset ) {
        _colorset = colorset;
        for( final TreePanel p : getTreePanels() ) {
            p.setBackground( colorset.getBackgroundColor() );
        }
    }

    void setTreeFontSet( final TreeFontSet fontset ) {
        _fontset = fontset;
    }

    void setupTreeGraphic( final Configuration config_settings, final ControlPanel control ) {
        control.setSpeciesColors( config_settings.getSpeciesColors() );
        control.setAnnotationColors( config_settings.getAnnotationColors() );
        RenderableDomainArchitecture.setColorMap( config_settings.getDomainColors() );
    }

    public void setWaitCursor() {
        setCursor( TreePanel.WAIT_CURSOR );
        repaint();
    }

    void terminate() {
        for( final TreePanel atvtreepanel : _treepanels ) {
            atvtreepanel.removeAllEditNodeJFrames();
        }
    }
}
