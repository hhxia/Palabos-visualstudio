#pragma once

class CLBMAStrategy
{
private:
	CLBMAStrategy(void);

public:
	virtual ~CLBMAStrategy(void);

	CLBMAStrategy * GetNewInstance();
};
