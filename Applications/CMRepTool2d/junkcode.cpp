
static const int TOOL_NONE;
static const int TOOL_ZOOM;
static const int TOOL_PAN;
static const int TOOL_SELECTION;
static const int TOOL_TRANSLATE;
static const int TOOL_ROTATE;
static const int TOOL_SCALE;
static const int TOOL_STRETCH;
static const int TOOL_BEND;


/**
 * A tool object.  The tool handles input for a given tool.  
 */
class ToolHandler {
	DisplayWindow &dw;
public:
	virtual void click(MouseEvent &ev) {};
	virtual void drag(MouseEvent &ev) {};
	virtual void display() {};

	ToolHandler(DisplayWindow &win) : dw(win) {
	}
};

/**
 * Zoom tool.  
 */
class ZoomToolHandler : public ToolHandler {
public:
	void drag(MouseEvent &ev);

	ZoomToolHandler(DisplayWindow &win) : ToolHandler(win) {
	}
};


/**
 * Pan tool
 */
class PanToolHandler : public ToolHandler {
public:
	void drag(MouseEvent &ev);

	PanToolHandler(DisplayWindow &win) : ToolHandler(win) {
	}
};


/**
 * Selection tool.  
 */
class SelectionToolHandler : public ToolHandler {
public:
	void click(MouseEvent &ev);
	void drag(MouseEvent &ev);
	void display();

	SelectionToolHandler(DisplayWindow &win) : ToolHandler(win) {
	}
};


/**
 * Translate tool.  
 */
class TranslateToolHandler : public ToolHandler {
public:
	void drag(MouseEvent &ev);

	TranslateToolHandler(DisplayWindow &win) : ToolHandler(win) {
	}
};


/**
 * Rotate tool.  
 */
class RotateToolHandler : public ToolHandler {
public:
	void click(MouseEvent &ev);
	void drag(MouseEvent &ev);
	void display();

	RotateToolHandler(DisplayWindow &win) : ToolHandler(win) {
	}
};


/**
 * Scale tool.  
 */
class ScaleToolHandler : public ToolHandler {
public:
	void click(MouseEvent &ev);
	void drag(MouseEvent &ev);
	void display();

	ScaleToolHandler(DisplayWindow &win) : ToolHandler(win) {
	}
};


/**
 * Stretch tool.  
 */
class StretchToolHandler : public ToolHandler {
public:
	void click(MouseEvent &ev);
	void drag(MouseEvent &ev);
	void display();

	StretchToolHandler(DisplayWindow &win) : ToolHandler(win) {
	}
};


/**
 * Zoom tool.  
 */
class ToolHandler : public ToolHandler {
public:
	void click(MouseEvent &ev);
	void drag(MouseEvent &ev);
	void display();

	ZoomToolHandler(DisplayWindow &win) : ToolHandler(win) {
	}
};


/**
 * Zoom tool.  
 */
class ZoomToolHandler : public ToolHandler {
public:
	void click(MouseEvent &ev);
	void drag(MouseEvent &ev);
	void display();

	ZoomToolHandler(DisplayWindow &win) : ToolHandler(win) {
	}
};


/**
 * Zoom tool.  
 */
class ZoomToolHandler : public ToolHandler {
public:
	void click(MouseEvent &ev);
	void drag(MouseEvent &ev);
	void display();

	ZoomToolHandler(DisplayWindow &win) : ToolHandler(win) {
	}
};

