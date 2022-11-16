import os, csv, sqlite3

conn = sqlite3.connect('dados.db')
cursor = conn.cursor()

cursor.execute(
    """
        CREATE TABLE IF NOT EXISTS resultados (
            id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
            tempos VARCHAR(11) NOT NULL,
            tratamentos_evidencia VARCHAR(11) NOT NULL,
            tratamentos_comparados VARCHAR(11) NOT NULL,
            value VARCHAR(100) NOT NULL,
            str_error VARCHAR(100) NOT NULL,
            df VARCHAR(11) NOT NULL,
            t_value VARCHAR(100) NOT NULL,
            p_value VARCHAR(100) NOT NULL,
            significancia VARCHAR(5) NOT NULL
        )
    """
)

pasta_resultado = f'{os.getcwd()}/result'
for i, j in enumerate(os.listdir(pasta_resultado)):
    # print(os.path.join(pasta_resultado, j))
    caminho_completo = os.path.join(pasta_resultado, j)
    with open(caminho_completo, 'r', encoding='utf-8') as arquivo:
        leitor = csv.reader(arquivo, delimiter=',')
        for y, coluna in enumerate(leitor):
            if y == 0:
                continue
            
            if len(coluna) < 9:
                continue
            
            if i > 0 :
                
                try:
                    cursor.execute(
                        f"""
                        SELECT * FROM resultados WHERE tempos='{coluna[0]}' AND tratamentos_evidencia='{coluna[1]}' AND tratamentos_comparados='{coluna[2]}'
                        """
                    )
                    for linha in cursor.fetchall():
                        if linha:
                            # print('***************')
                            continue

                    cursor.execute(
                        f"""
                        SELECT * FROM resultados WHERE tempos='{coluna[0]}' AND tratamentos_evidencia='{coluna[2]}' AND tratamentos_comparados='{coluna[1]}'
                        """
                    )
                    for linha in cursor.fetchall():
                        if linha:
                            # print('***************')
                            continue
                    
                    cursor.execute(
                        f"""
                        INSERT INTO resultados (tempos, tratamentos_evidencia, tratamentos_comparados, value, str_error, df, t_value, p_value, significancia)
                        VALUES
                        ('{coluna[0]}', '{coluna[1]}', '{coluna[2]}', '{coluna[3]}', '{coluna[4]}', '{coluna[5]}', '{coluna[6]}', '{coluna[7]}', '{coluna[8]}')
                        """
                    )
                except Exception as ex:
                    print(ex) 
    
conn.commit()


try:
    cursor.execute(
        f"""
        SELECT * FROM resultados
        """
    )
    for linha in cursor.fetchall():
        # id, tempos, tratamentos_evidencia, tratamentos_comparados, value, str_error, df, t_value, p_value, significancia = linha
        print(linha)
except Exception as ex:
    print(ex)

# input('Marlon: ')

try:
    cursor.execute(
        f"""
        SELECT * FROM resultados
        """
    )
    f = open(os.path.join(pasta_resultado, 'resultados.csv'),'w', newline='', encoding='utf-8')
    w = csv.writer(f)
    w.writerow(['ID', 'TEMPOS', 'TRATAMENTO EM EVIDENCIA', 'TRATAMENTO COMPARADO', 'VALUE', 'STD.ERROR', 'DF', 'T.VALUE', 'P.VALUE', 'SIGN.'])
    for linha in cursor.fetchall():
        # id, tempos, tratamentos_evidencia, tratamentos_comparados, value, str_error, df, t_value, p_value, significancia = linha
        if len(linha) == 10:
            # print(linha)
            w.writerow(linha)
except Exception as ex:
    print(ex)

os.remove(os.path.join(os.getcwd(), 'dados.db'))