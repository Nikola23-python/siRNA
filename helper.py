class helper:
    """Документация рабочая"""
    def info_about_methods(obj):
        """Функция для получения информации о методах класса"""
        print("=" * 60)
        print(f"ДОКУМЕНТАЦИЯ КЛАССА: {obj.__class__.__name__}")
        print("=" * 60)
        class_doc = obj.__class__.__doc__
        if class_doc:
            print(f"📚 Документация класса:\n{class_doc.strip()}\n")
        else:
            print("📚 Документация класса: Отсутствует\n")
        methods = [method for method in dir(obj)
                   if not method.startswith('_') or
                   (method.startswith('_') and not method.startswith('__'))]
        for method_name in methods:
            try:
                method = getattr(obj, method_name)
                if callable(method):
                    print(f"🔹 МЕТОД: {method_name}()")
                    doc = method.__doc__
                    if doc:
                        lines = doc.strip().split('\n')
                        for line in lines:
                            print(f"   {line.strip()}")
                    else:
                        print(" 📝 Документация: Отсутствует")
                    print()
            except Exception as e:
                print(f"❌ Ошибка при обработке метода {method_name}: {e}")



    """Основной метод проекта, но пока без понятия как реализовать(не работает)"""
    def run_all_methods(self):
        """Автоматически находит и запускает все необходимые методы класса"""
        exclude_methods = ['info_about_methods']
        results = []
        methods = [method for method in dir(self)
                   if not method.startswith('_') and
                   callable(getattr(self, method)) and
                   method != 'run_all_methods' and
                   method not in exclude_methods]

        print(f"📋 Найдено методов: {len(methods)}")

        for method_name in methods:
            try:
                method = getattr(self, method_name)
                print(f"▶️  Выполняется: {method_name}()")
                result = method()
                results.append({
                    'method': method_name,
                    'result': result,
                    'doc': method.__doc__.strip() if method.__doc__ else 'Без описания'
                })
                print(f"   ✅ Результат: {result}")

            except Exception as e:
                results.append({
                    'method': method_name,
                    'result': f"Ошибка: {e}",
                    'doc': 'Метод завершился с ошибкой'
                })
                print(f" ❌ Ошибка: {e}")
        return results

    def to_csv(self, filename=None):
        if filename:
            self.to_csv(filename, index=False)
            print(f"Результаты сохранены в файл: {filename}")
        return self